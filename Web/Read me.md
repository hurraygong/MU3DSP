1. Introduction

MU3DSP, a starting from protein sequence prediction tool to assess single point mutation effects on protein thermodynamic stability and applying to dingle-domain monomeric proteins. Given protein sequence with single mutations as the input, the proposed model integrated both sequence level features of mutant residues and residue level mutation-based 3D structure features. 

**Figure 1**

2. Inputs

MU3DSP provides two predicting modes, normal mode and Fast mode. Normal mode can achieve more precision result than fast mode, but need more time. The normal mode is recommended. However, if you don't want to using structures of the query variants or even their homologous protein structures, the fast mode will be best.


MU3DSP supports processing two input data format:

a). protein sequence 

The input format is:
```
SEQ: SEQID Sequence
MUT: <wild-type1><position1><mutation1> <wild-type2><position2><mutation2> <wild-type3><position3><mutation3>
```

- Input files should not contain suffixes.
- Each mutation to be predicted contains two lines, the first line should start with 'SEQ:' and the second line should start with 'MUT:'.
- 'SEQ:' line defines the target protein sequence in a one-letter code: 'SEQ: SEQId Sequence'. Each element is split by space.
- 'MUT:' line defines the mutations which are to be predicted for the protein sequence specified on the previous line: 'MUT: <wild-type><position><mutation> ...'. Each mutation is split by space.
- proteinId: protein sequence unique identifier (the identifier should not contain '.'), It is better to use Uniprot ID or PDB ID added chain ID.
- mutation: The mutation format '<wild-type><position><mutation>' stands for wild-type amino acid, mutation position in sequence, and mutated amino acid.
- Multiple mutations per line are allowed, however, they are always predicted as single-site mutations at a time.
- Please note that unknown amino acids such as 'X' should be removed beforehand.
- The whole protein sequence must be specified on a single line, i.e., line-breaks are not allowed.

For example:
```
SEQ: 2ocj_A SVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENL
MUT: Q9H Q9P

SEQ: 2ocj_A SVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENL
MUT: Q9*
```
    
b). protein sequence with `.a3m` file from `HHblits`.

we also provide the interface `HHblits`, the upload file must be MSA file of the input sequece

you can input your Email address for receiving results.

Running time is related to the homologous protein structures, the running time of MSA and the condition of our server. For example, if the structures are unavailable and the MSA file(.a3m format) are provided, it will take less 1 minute to finish. Therefore, it better to input one file for mutations in one protein. 


3. Outputs

We will send the results to your email when the job is finished. Results will be shown in the result page (example) when the job is finished. In addition, results can be downloaded by clicking "Download results". 