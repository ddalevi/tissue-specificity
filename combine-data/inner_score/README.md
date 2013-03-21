In order to calculate inner score, run the script inner_score.pl on the output files of each data set. 
The outcomes can then be used to calculate the combined total score. 

linux> ./inner_score.pl *_result.txt *_uniq.txt outputfile.txt

1. inner_score.pl combines probesets and map to the gene symbols, then calculate inner score in one data set.
2. *_result.txt is the result from each data set *. 
	serial_nr	probeset	gene_symbol	nr_specific_tissue	tissue1	tissue2
	10	1431_at	CYP2E1	1	Liver	NA
3. *_uniq.txt is a list of non-duplicated gene symbols in data set *.
4. specify your outputfile so the outcome will be stored in outputfile.txt


E.g. calculate the inner score of GSE7307, 
linux> ./inner_score.pl gse7307_result.txt gse7307_uniq.txt inner_gse7307.txt

The outcome file inner_gse7307.txt lists:
	gene_symbol	tissue_name1	tissue_score1	tissue_name2(if any)	tissue_score2(if any)	END
e.g.	USP6NL			Esophagus	1	END
	CNGA3			SpinalCord	0.5	CNS	0.5	END

