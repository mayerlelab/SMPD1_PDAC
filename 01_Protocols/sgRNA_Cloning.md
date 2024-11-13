# Cloning sgRNA into Vector

## Introduction
This protocol outlines the step-by-step process for inserting sgRNA into the CRISPR/Cas9 vector pSpCas9(BB)-2A-Puro (PX459) V2.0. The sgRNAs are designed to target specific regions within a gene of interest.

## Materials

### Plasmid
- pSpCas9(BB)-2A-Puro (PX459) V2.0

### sgRNA Oligos
- Custom-designed to target specific gene regions
    - Forward: 5’ <r>CACC</r><g>G</g><bl>NNNNNNNNNNNNNNNNNNNN</bl>3’

        <r>CACC: overhang for BbsI cloning site </r>

        <g>G: for start at U6 expression site (only if gRNA start does not have G) </g>

        <bl>N: gRNA sequence (forward) </bl>

    - Reverse: 5’ <o>AAAC</o><r>RRRRRRRRRRRRRRRRRRRR</r><g>C</g>3’

        <o>AAAC: overhang for BbsI cloning site</o>

        <g>C: for start at U6 expression site (only if gRNA start does not have G)</g>

        <r>R: gRNA sequence (Reverse complement)<r>

### Buffers
- Tango buffer (10x)
- T4 ligation buffer (10x)
- DTT (10 mM)
- ATP (10 mM)
- Nuclease-free water

### Enzymes
- <i>Bbsi</i>
- T4 Polynucleotide Kinase (PNK)
- T7 Ligase

## Procedure

### Preparation of the sgRNA Oligos Inserts

1. Resuspend the sgRNA oligos to a final concentration of 100 µM.

1. phosphorylation and annealing: 

| **Component**           | **Amount (µl)** |
|-------------------------|-----------------|
| sgRNA Forward (100 µM)  | 1               |
| sgRNA Reverse (100 µM)  | 1               |
| T4 ligation buffer 10x  | 1               |
| T4 PNK                  | 1               |
| ddH2O                   | 6               |
| **Total**               | **10**          |

1. Phosphorylate and anneal the oligos:  

   Use a thermocycler with the following program:  
   - 37°C for 30 minutes  
   - 95°C for 5 minutes  
   - Ramp down to 25°C at 5°C/min

1. Dilute the phosphorylated and annealed oligos:  
   - Dilute the oligos 1:200 in nuclease-free water at room temperature.

### Ligation of the Phosphorylated and Annealed sgRNA Oligos into pSpCas9(BB)-2A-Puro (PX459) V2.0


1. Prepare the following ligation mixture:  

| **Component**            | **Amount (µl)** |
|--------------------------|-----------------|
| PX459, 100 ng            | x               |
| Diluted oligo duplex     | 2               |
| Tango buffer, 10x        | 2               |
| DTT, 10 mM               | 1               |
| ATP, 10 mM               | 1               |
| FastDigest <i>BbsI</i>   | 1               |
| T7 ligase                | 0.5             |
| ddH2O                    | To 20           |

1. Incubate the ligation reaction:  
   - Incubate the mixture in PCR  with following cycling condtions.

| **Cycle number** | **Condition**                 |
|------------------|-------------------------------|
| 1-6              | 37°C for 5 min, 21°C for 5 min |

### Transformation

1. Thaw DH5α cells on ice for 30 minutes.
1. Add 5 µL of the ligation reaction to 70 µL of DH5α competent cells in a new 1.5 mL tube.
1. Incubate the tube on ice for 12 minutes.
1. Heat shock the cells at 42°C for 2 minutes.
1. Immediately place the tube back on ice for 2 minutes.
1. Add 100 µL of LB medium to the tube and incubate at 37°C with shaking for 45 minutes.
1. Plate the aliquots on LB agar plates containing the appropriate antibiotic.
1. Incubate the plate overnight at 37°C.

### Mini Prep

1. Following morning, pick 3 different colonies from each plate.

1. Incubate each colony in 8 mL of LB medium containing 100 µg/mL ampicillin in 10 mL polyethylene tubes 6 hours at 37°C with shaking at 150 rpm

1. Following incubation, follow the procedure outlined in the HiYield Plasmid Mini Kit to extract the plasmid DNA.

1. Measure the DNA concentration using a spectrophotometer.

1. Prepare 100 ng of plasmid DNA in a total volume of 15 µL nuclease-free water.

1. Add 2 µL of 10 µM sequencing primer to the plasmid DNA.

1. Label the samples appropriately and send the sequencing reactions to Eurofins for analysis.
