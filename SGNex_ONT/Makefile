all_preds: clean
	./organize_fusion_preds.pl
	./SGNEx-analyze_Illum_as_TP_unique_as_FP.pl


restrict_fuzzy_breakpoints: clean
	./organize_fusion_preds.pl --restrict_fuzzy_breakpoints
	./SGNEx-analyze_Illum_as_TP_unique_as_FP.pl


clean:
	./cleanMe.sh
