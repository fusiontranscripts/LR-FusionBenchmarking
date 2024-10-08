all_preds: clean
	./organize_fusion_preds.pl
	./SGNEx-analyze_Illum_as_TP_unique_as_FP.pl preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments either
	./SGNEx-analyze_Illum_as_TP_unique_as_FP.pl preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments starF
	./SGNEx-analyze_Illum_as_TP_unique_as_FP.pl preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments arriba
	./SGNEx-analyze_Illum_as_TP_unique_as_FP.pl preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments starF,arriba

restrict_fuzzy_breakpoints: clean
	./organize_fusion_preds.pl --restrict_fuzzy_breakpoints
	./SGNEx-analyze_Illum_as_TP_unique_as_FP.pl preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.FUZZY_RESTRICTED either
	./SGNEx-analyze_Illum_as_TP_unique_as_FP.pl preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.FUZZY_RESTRICTED starF
	./SGNEx-analyze_Illum_as_TP_unique_as_FP.pl preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.FUZZY_RESTRICTED arriba
	./SGNEx-analyze_Illum_as_TP_unique_as_FP.pl preds.collected.gencode_mapped.wAnnot.filt.pass.proxy_assignments.FUZZY_RESTRICTED starF,arriba

clean:
	./cleanMe.sh

