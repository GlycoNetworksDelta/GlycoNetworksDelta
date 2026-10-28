Glycosylation is critical to many biological processes, such as protein folding and intercellular communication. Mass spectrometry is widely used for studying proteins and glycoproteins. For glycoproteins, successful data interpretation of mass spectrometry data using database search requires a complete glycan database — a condition that is often not satisfied.

GlycoNetworks builds improved N-glycan databases by constructing a sample-specific glycan database based on the mass spectrometry data itself, rather than relying solely on existing glycan databases. GlycoNetworks reads in mass spectrometry data and an initial (seed) glycan database. Then it analyzes the data to augment the initial glycan database with additional glycans. Rather than analyzing individual spectra in isolation, GlycoNetworks constructs a network (graph) of putative glycans where each node represents a glycan and each edge represents a monosaccharide. More details can be found in the paper ["Network method for building a sample-specific glycan database for N-linked glycosylation from MS/MS data."](https://doi.org/10.1021/acs.jproteome.5c00006)

A limitation of GlycoNetworks is that the user must specify all building blocks (such as monosaccharides or modifications) in advance. GlycoNetworksΔ extends GlycoNetworks and removes this limitation, enabling the discovery of unanticipated modifications to N-glycans. More details can be found in the paper "Finding Unanticipated Modifications to N-Glycans in Glycoproteomics Data."

## Usage summary

GlycoNetworksΔ is a command line tool.
```
Usage: GlycoNetworksDelta <Spectra file path> <Glycan database file path> [Mass tolerance] [Deltas]
  <Spectra file path> is mandatory. Must be in MGF format, or ConvertToMgf.bat must be configured to convert to MGF format
  <Glycan database file path> is mandatory. Must be a text file of glycan compositions
  [Verbosity] is optional. Current verbosity levels are 0, 1, or 2. Level 0 outputs delta mass and confidence interval. Level 1 provides basic detail. Level 2 shows 'all' details. Default is: 1
  [Mass tolerance] is optional. Default is: 0.03
  [Deltas] is optional. List of glycan deltas to consider, separated by commas. You may also specify a delta mass
    Default is: all building blocks found in the glycan database, along with the disaccharide Hex-HexNAc
    This program knows the following deltas: Hex,HexNAc,Fuc,NeuAc,NeuGc,Pent,Phospho,Na,Acetyl,Sulfo,pC,Hex-HexNAc
  ```

## Usage example

The `GlycanTestData` folder has example data for testing. This data is from dataset [PXD041208](https://www.ebi.ac.uk/pride/archive/projects/PXD041208) in the PRIDE data repository and is discussed in [this paper](https://doi.org/10.1038/S41467-024-46443-0). To reduce file size, only a narrow retention-time window, in which a significant proportion of glycopeptides elute, has been extracted from the C1 cell line antibody sample acquired on the Orbitrap Exploris 480.

Example command for running GlycoNetworks:
```
GlycoNetworksDelta.exe 20220920_EX3_UM7_shamo002_SA_EXT02_C1_Abs_pd.12-16min.mgf "N-glycan 59 common biantennary.txt" 1 0.03 HexNAc,Fuc,Hex,NeuAc,NeuGc,Hex-HexNAc
```

Console output:
```
Deltas: HexNAc,Fuc,Hex,NeuAc,NeuGc,Hex-HexNAc

Number of spectra: 2822
Number of spectra with glycopeptide oxonium ion and Y0, Y1, Y2, etc.: 499 (17.68%)

Additional glycans:
HexNAc(2)Hex(10)
HexNAc(2)Hex(11)
HexNAc(3)Hex(10)
HexNAc(3)Hex(11)
HexNAc(3)Hex(12)
HexNAc(3)Hex(7)
HexNAc(3)Hex(8)
HexNAc(3)Hex(9)

Candidate delta masses:
Mass: 27.9952   CI: (27.9943, 27.9960)          Quality score: 48.3
Mass: 63.9157   CI: (63.9146, 63.9167)          Quality score: 14.1
```

The top candidate delta mass is likely formylation (probably a sample prep artifact, as discussed in [this paper](https://doi.org/10.1007/s00216-022-03989-6)). We can run GlycoNetworksΔ again to construct a glycan database with the formylation:
```
GlycoNetworksDelta.exe 20220920_EX3_UM7_shamo002_SA_EXT02_C1_Abs_pd.12-16min.mgf "N-glycan 59 common biantennary.txt" 1 0.03 HexNAc,Fuc,Hex,NeuAc,NeuGc,Hex-HexNAc,27.994915
```

Console output:
```
Deltas: HexNAc,Fuc,Hex,NeuAc,NeuGc,Hex-HexNAc,27.994915

Number of spectra: 2822
Number of spectra with glycopeptide oxonium ion and Y0, Y1, Y2, etc.: 499 (17.68%)

Additional glycans:
HexNAc(2)Hex(10)
HexNAc(2)Hex(11)
HexNAc(2)Hex(2) 27.994915
HexNAc(2)Hex(3) 27.994915
HexNAc(2)Hex(4) 27.994915
HexNAc(2)Hex(5) 27.994915
HexNAc(2)Hex(6) 27.994915
HexNAc(2)Hex(7) 27.994915
HexNAc(2)Hex(8) 27.994915
HexNAc(2)Hex(9) 27.994915
HexNAc(3)Hex(10)
HexNAc(3)Hex(11)
HexNAc(3)Hex(12)
HexNAc(3)Hex(3) 27.994915
HexNAc(3)Hex(4) 27.994915
HexNAc(3)Hex(5) 27.994915
HexNAc(3)Hex(6) 27.994915
HexNAc(3)Hex(7)
HexNAc(3)Hex(8)
HexNAc(3)Hex(9)

Candidate delta masses:
Mass: 63.9157   CI: (63.9146, 63.9168)          Quality score: 17.8
Mass: 1159.5058 CI: (1159.5048, 1159.5068)      Quality score: 3.7
```
