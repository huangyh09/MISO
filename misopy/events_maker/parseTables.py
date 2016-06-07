import os, sys, operator, string
import collections

# import rnaseqlib
# import rnaseqlib.tables as tables

# def table_fname_to_header(table_fname):
#     """
#     Returns the header for the table.
#     """
#     name = os.path.basename(table_fname)
#     header = None
#     if name.startswith("ensGene"):
#         header = tables.UCSC_ENSGENE_HEADER
#     elif name.startswith("refGene"):
#         header = tables.UCSC_REFGENE_HEADER
#     elif name.startswith("knownGene"):
#         header = tables.UCSC_KNOWNGENE_HEADER
#     else:
#         return None
#     # Designate one of the columns as a 'gene' column. Used to
#     # annotate which event a gene falls in.
#     if "name2" in header:
#         # If name2 exists (like it does in refGene and ensGene
#         # tables), use it as the gene column
#         header[header.index("name2")] = "gene"
#     elif "gene" in header:
#         pass
#     else:
#         # If name2 isn't present (e.g. in knownGene), then
#         # use the transcript ID field 'name'
#         header[header.index("name")] = "gene"
#     return header
    

# # Generic function to read in a file.
# def readTable(table_f):
#     data = [] 
#     colToIdx = {}

#     table_in = open(table_f)
#     header = table_fname_to_header(table_f)
#     if header is None:
#         raise Exception, "Unrecognized table file %s" %(table_f)

#     for line in open(table_f):
#         if line.startswith("#"):
#             # Skip header
#             continue
#         else:
#             vals = line.strip().split("\t")
#             # Mapping from column names in UCSC table to values
#             # when no header is given
#             col_values = \
#                 dict([(header[col_num], vals[col_num]) \
#                       for col_num in range(len(header))])
#             item = [col_values["chrom"],
#                     col_values["exonStarts"],
#                     col_values["exonEnds"],
#                     col_values["strand"],
#                     col_values["gene"]]
#             data.append(item)
#     return data



# Generic function to read in a file in UCSC table format.
def readTable_ucsc(table_f):
    """Make sure that the UCSC table has chrom, strand, exonStarts, 
    exonEnds, and geneID in the right columns.
    """
    data = []
    for line in open(table_f):
        if line.startswith("#"):
            # Skip header
            continue
        else:
            vals = line.strip().split("\t")
            item = [vals[2], vals[9], vals[10], vals[3], vals[12]]
            data.append(item)
    return data


# Generic function to read in a file in gff format.
# This change included by Yuanhua Huang <Y.Huang@ed.ac.uk>
def readTable_gff(table_f, ftype="gff3"):
    """We assume that the gene comes out before transcript;
    and transcript comes out before exons.
    """
    fid = open(table_f, "r")
    all_lines = fid.readlines()
    fid.close()

    data = []
    for i in range(len(all_lines)):
        line = all_lines[i]
        if line.startswith("#") or line.startswith(">"): 
            # Skip header
            continue

        vals = line.strip().split("\t")
        if vals[2] == "gene":
            gene_id = "*"
            if ftype == "gtf":
                idx = vals[8].find("gene_id")
                if idx>-1: gene_id = vals[8][idx:].split('"')[1]
            if ftype == "gff3":
                idx = vals[8].find("ID")
                if idx>-1: gene_id = vals[8][idx+3:].split(';')[0]
            continue

        elif vals[2] == "transcript" or vals[2] == "mRNA":
            chrom = vals[0]
            strand = vals[6]
            exonStarts = []
            exonStops = []
            continue

        elif vals[2] == "exon":
            exonStarts.append(str(int(vals[3])-1))
            exonStops.append(vals[4])
            
        if (i == len(all_lines)-1 or 
            all_lines[i+1].strip().split("\t")[2] == "gene" or
            all_lines[i+1].strip().split("\t")[2] == "transcript" or
            all_lines[i+1].strip().split("\t")[2] == "mRNA"):

            exonStarts = ",".join(exonStarts)
            exonStops = ",".join(exonStops)
            data.append([chrom, exonStarts, exonStops, strand, gene_id])
    return data


def readTable(table_f, ftype="gtf"):
    """read table or gtf/gff3 files.
    """
    if ftype == "gff3" or ftype == "gtf":
        return readTable_gff(table_f, ftype)
    elif ftype == "ucsc":
        return readTable_ucsc(table_f)


# Get splice graph.
def populateSplicegraph(table_f, ss5_ss3_F, ss3_ss5_F, ss5_ss3_R, ss3_ss5_R):
    data = readTable(table_f)
  
    #ss5_ss3_F = {}    # donor to acceptor (forward)
    #ss3_ss5_F = {}    # acceptor to donor (forward)
    #ss5_ss3_R = {}    # donor to acceptor (reverse)
    #ss3_ss5_R = {}    # acceptor to donor (reverse)
 
    for item in data: 
        chromval, startvals, endvals, strandval, gene = item
        startvals = map(int, startvals.split(",")[:-1])
        # Adds +1 since downloaded UCSC tables are 0-based start!
        startvals = map(str, [x + 1 for x in startvals])
        endvals = endvals.split(",")[:-1]
        if strandval == '+':
            for i in range(len(startvals) - 1):
                prevacceptor = ":".join([chromval, startvals[i], strandval])
                donor = ":".join([chromval, endvals[i], strandval])
                acceptor = ":".join([chromval, startvals[i + 1], strandval])
                nextdonor = ":".join([chromval, endvals[i + 1], strandval])

                if donor not in ss5_ss3_F:
                    ss5_ss3_F[donor] = []
                ss5_ss3_F[donor].append(acceptor)
                if acceptor not in ss3_ss5_F:
                    ss3_ss5_F[acceptor] = []
                ss3_ss5_F[acceptor].append(nextdonor)

                if donor not in ss5_ss3_R:
                    ss5_ss3_R[donor] = []
                ss5_ss3_R[donor].append(prevacceptor)
                if acceptor not in ss3_ss5_R:
                    ss3_ss5_R[acceptor] = []
                ss3_ss5_R[acceptor].append(donor)

        else:
            startvals = startvals[::-1]
            endvals = endvals[::-1]
            for i in range(len(startvals) - 1):
                prevacceptor = ":".join([chromval, endvals[i], strandval])
                donor = ":".join([chromval, startvals[i], strandval])
                acceptor = ":".join([chromval, endvals[i + 1], strandval])
                nextdonor = ":".join([chromval, startvals[i + 1], strandval])

                if donor not in ss5_ss3_F:
                    ss5_ss3_F[donor] = [] 
                ss5_ss3_F[donor].append(acceptor)
                if acceptor not in ss3_ss5_F:
                    ss3_ss5_F[acceptor] = []
                ss3_ss5_F[acceptor].append(nextdonor)

                if donor not in ss5_ss3_R:
                    ss5_ss3_R[donor] = []
                ss5_ss3_R[donor].append(prevacceptor)
                if acceptor not in ss3_ss5_R:
                    ss3_ss5_R[acceptor] = []
                ss3_ss5_R[acceptor].append(donor)

    return ss5_ss3_F, ss3_ss5_F, ss5_ss3_R, ss3_ss5_R 


# Get counts of each splice site.    
def cleanSplicegraph(ss5_ss3_F, ss3_ss5_F, ss5_ss3_R, ss3_ss5_R):

    for ss in ss5_ss3_F:
        ss5_ss3_F[ss] = collections.Counter(ss5_ss3_F[ss])
    for ss in ss3_ss5_F:
        ss3_ss5_F[ss] = collections.Counter(ss3_ss5_F[ss])
    for ss in ss5_ss3_R:
        ss5_ss3_R[ss] = collections.Counter(ss5_ss3_R[ss])
    for ss in ss3_ss5_R:
        ss3_ss5_R[ss] = collections.Counter(ss3_ss5_R[ss])

    return ss5_ss3_F, ss3_ss5_F, ss5_ss3_R, ss3_ss5_R 


def readXref(xref_f):
  
    geneToInfo = {} 
    for line in open(xref_f): 
        
        tx, gene, symbol, desc = line.strip().split("\t")
        if gene not in geneToInfo:
            geneToInfo[gene] = [symbol, desc]
        else:
            if symbol != 'n/a' and geneToInfo[gene][0] == 'n/a':
                geneToInfo[gene][0] = symbol
            if desc != 'n/a' and geneToInfo[gene][1] == 'n/a':
                geneToInfo[gene][1] = desc
    return geneToInfo


def populateGenelist(table_f):
   
    data = readTable(table_f)
    ssToGene = {}

    for item in data: 
        chromval, startvals, endvals, strandval, gene = item
        startvals = map(int, startvals.split(",")[:-1])
        startvals = map(str, [x + 1 for x in startvals])
        endvals = endvals.split(",")[:-1]
        for i in range(len(startvals)):
            ss1 = ":".join([chromval, startvals[i], strandval])
            ss2 = ":".join([chromval, endvals[i], strandval])
            ssToGene[ss1] = gene
            ssToGene[ss2] = gene
    
    return ssToGene 




