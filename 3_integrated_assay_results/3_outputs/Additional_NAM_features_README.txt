plddt_f32: AlphaFold per-residue confidence score (0-100). Low pLDDT often indicates intrinsically disordered regions (IDRs). 
		Variants in IDRs behave differently than in structured domains.

rsa_f32: Relative Solvent Accessibility. Distinguishes buried residues (hydrophobic core, critical for stability) from surface residues 
		(interaction interfaces). Ranges from 0 (buried) to 100 (exposed)

AlphaMissense_score: AlphaMissense is an unsupervised model for predicting the pathogenicity of human missense variants
		by incorporating the structural context of an AlphaFold-derived system. The AlphaMissense score ranges from 0 to 1. 
		The larger the score, the more likely the variant is pathogenic. Details see https://doi.org/10.1126/science.adg7492.

ESM1b_score: ESM1b scores are log-likelihood ratio (LLR) scores for predicting the pathogenic effects of coding
		variants based on a 650-million-parameter protein language model, ESM1b. The range of ESM1b score in dbNSFP 
		is from -24.538 to 6.937. The smaller the score, the more likely the variant is pathogenic. 
		Details see doi: https://doi.org/10.1038/s41588-023-01465-0

REVEL_score: REVEL is an ensemble score based on 13 individual scores for predicting the
		pathogenicity of missense variants. Scores range from 0 to 1. The larger the score the more 
		likely the SNP has a damaging effect. "REVEL scores are freely available for non-commercial use.  
		For other uses, please contact Weiva Sieh" (weiva.sieh@mssm.edu)

MutPred2_score: General MutPred2 score. Scores range from 0 to 1. The larger the score the more
		likely the SNP has a damaging effect.

gnomAD4.1_joint_nhomalt: Count of individuals with homozygous alternative allele in the whole gnomAD joint (genome+exome) samples v4.1

phyloP100way_vertebrate: phyloP (phylogenetic p-values) conservation score based on the
		multiple alignments of 100 vertebrate genomes (including human). The larger the score, 
		the more conserved the site. Scores range from -20.0 to 10.003 in dbNSFP.

phastCons100way_vertebrate: phastCons conservation score based on the multiple alignments
		of 100 vertebrate genomes (including human). The larger the score, the more conserved 
		the site. Scores range from 0 to 1. 

GERP++_RS: GERP++ RS score, the larger the score, the more conserved the site. Scores range from
		-12.3 to 6.17.

PROVEAN_score: PROVEAN score (PROVEANori). Scores range from -14 to 14. The smaller the score
		the more likely the SNP has damaging effect. 
  
DSSP provides an elaborate description of the secondary structure elements in a protein structure, 
		including backbone hydrogen bonding and the topology of β-sheets. 
		The most popular feature is the per-residue assignment of secondary structure with a single character code:

		H = α-helix
		B = residue in isolated β-bridge
		E = extended strand, participates in β ladder
		G = 310-helix
		I = π-helix
		P = κ-helix (poly-proline II helix)
		T = hydrogen-bonded turn
		S = bend
		C = Coil

		g2p_Secondary structure (DSSP 9-state)* lists all 9 possible state 
		while g2p_Secondary structure (DSSP 3-state)* only lists 3, H (Helix), B (Beta sheet), C (Coil)
