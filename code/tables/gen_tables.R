#!/usr/bin/env Rscript
## Tables 1 (ACMG cutoffs, from model/) and S7 (gnomAD distribution, shipped); S4/S5
## (benchmark + leakage) need the HGMD-licensed test set under dataset/test_data/, else skipped.
suppressPackageStartupMessages({library(data.table); library(pROC); library(openxlsx)})
ROOT  <- Sys.getenv("INDELVAR_ROOT",".")
MODEL <- file.path(ROOT,"model")
HELD  <- file.path(ROOT,"dataset","test_data")   # HGMD-licensed test set (user-supplied)
PAPER <- file.path(ROOT,"dataset","table_fig_data")
TAB   <- file.path(ROOT,"tables"); dir.create(TAB, showWarnings=FALSE)
say   <- function(...) cat(sprintf(...),"\n")
wr <- function(dt, name){ fwrite(dt, file.path(TAB,paste0(name,".csv")))
  openxlsx::write.xlsx(dt, file.path(TAB,paste0(name,".xlsx"))); say("  wrote %s .csv/.xlsx (%d rows)", name, nrow(dt)) }

## ---------- Table 1: per-type cutoffs + evidence LR targets ----------
c1<-fread(file.path(MODEL,"evidence_cutoffs.tsv"))       # per-type ACMG cutoffs
ct<-fread(file.path(MODEL,"cutoffs_local_1sided95.tsv")) # evidence LR targets
ev<-c("PP3_Strong","PP3_Moderate3","PP3_Moderate","PP3_Supporting","BP4_Supporting","BP4_Moderate","BP4_Moderate3","BP4_Strong")
lab<-c("PP3 Strong (+4)","PP3 Moderate (+3)","PP3 Moderate (+2)","PP3 Supporting (+1)","BP4 Supporting (−1)","BP4 Moderate (−2)","BP4 Moderate (−3)","BP4 Strong (−4)")
gcut<-function(sub,e){v<-c1[subset==sub&evidence==e]$cutoff;if(length(v)==0||is.na(v))NA else v}
gtar<-function(sub,e){v<-ct[subset==sub&evidence==e]$target;if(length(v)==0||is.na(v))NA else v}
cut_str<-function(sub,e){v<-gcut(sub,e);if(is.na(v))"not reached" else sprintf("%s %.3f",ifelse(grepl("PP3",e),"≥","≤"),v)}
T1<-data.table(`ACMG/AMP stratum`=lab,
  `Deletion target LR`=sapply(ev,function(e)round(gtar("del",e),2)),
  `Deletion score cutoff`=sapply(ev,function(e)cut_str("del",e)),
  `Insertion target LR`=sapply(ev,function(e)round(gtar("ins",e),2)),
  `Insertion score cutoff`=sapply(ev,function(e)cut_str("ins",e)))
wr(T1,"table_1_acmg_cutoffs")

## ---------- Tables S4 / S5: test-set benchmark + leakage sensitivity (HGMD-licensed; skipped if absent) ----------
S4F<-file.path(HELD,"test_set_tool_benchmark_scored.tsv")
if(!file.exists(S4F)){
  say("[S4/S5] test-set benchmark not shipped (HGMD-licensed); build it with benchmark.R under your own HGMD license and place it under dataset/test_data/.")
} else {
  m<-fread(S4F); m[,y:=as.integer(label=="Pathogenic")]
  TOOLS<-c("INDELVAR","IndeLLM","FATHMM-indel","CADD","ESM-1b","ProGen2","phyloP100","GERP++","SHINE","PON-Del","MutPred-Indel")
  disp<-c(INDELVAR="INDELVAR",IndeLLM="IndeLLM *",`FATHMM-indel`="FATHMM-indel",CADD="CADD v1.7",`ESM-1b`="ESM-1b",
          ProGen2="ProGen2",phyloP100="phyloP100",`GERP++`="GERP++",SHINE="SHINE *",`PON-Del`="PON-Del *",`MutPred-Indel`="MutPred-Indel *")
  roc1<-function(y,s){r<-roc(y,s,quiet=T,levels=c(0,1),direction="<");c<-as.numeric(ci.auc(r));list(auc=as.numeric(auc(r)),lo=c[1],hi=c[3],roc=r)}
  orient<-function(s,y){ if(mean(s[y==1])<mean(s[y==0])) -s else s }
  rows<-list()
  for(ty in c("deletion","insertion")){d<-m[indel_type==ty];nc<-nrow(d)
    for(t in TOOLS){ if(!t%in%names(d))next
      xo<-d[!is.na(suppressWarnings(as.numeric(get(t))))]
      if(nrow(xo)<20||length(unique(xo$y))<2)next
      rr<-roc1(xo$y, orient(as.numeric(xo[[t]]),xo$y))            # AUROC on tool's own coverage
      if(t=="INDELVAR"){ p<-NA } else {                          # DeLong on INDELVAR-paired common subset
        xc<-d[!is.na(INDELVAR)&!is.na(suppressWarnings(as.numeric(get(t))))]
        p<-tryCatch(roc.test(
          roc(xc$y,xc$INDELVAR,quiet=T,levels=c(0,1),direction="<"),
          roc(xc$y,orient(as.numeric(xc[[t]]),xc$y),quiet=T,levels=c(0,1),direction="<"),
          paired=TRUE,method="delong")$p.value, error=function(e)NA) }
      rows[[length(rows)+1]]<-data.table(`Indel class`=sprintf("%s (n=%s)",ifelse(ty=="deletion","Deletions","Insertions"),format(nc,big.mark=",")),
        Tool=disp[t],AUROC=round(rr$auc,3),`95% CI`=sprintf("%.3f-%.3f",rr$lo,rr$hi),
        `Coverage %`=round(100*nrow(xo)/nc,1),`n scored`=nrow(xo),`n P`=sum(xo$y),`n B`=sum(1-xo$y),
        `DeLong p (vs INDELVAR)`=if(is.na(p))"" else signif(p,2))}}
  S5bench<-rbindlist(rows)
  ## Leakage-free DeLong (INDELVAR vs asterisk-marked tool) after removing each tool's
  ## training overlap (id lists under leakage/).
  LK<-file.path(HELD,"leakage")
  lfrow<-function(tool,cls,fmt="%.1f"){
    d<-if(cls=="all") m else m[indel_type==cls]
    ov<-readLines(file.path(LK,sprintf("%s_train_overlap_ids.txt",tolower(gsub("[^A-Za-z]","",tool)))))
    ov<-ov[nzchar(ov)]
    si<-suppressWarnings(as.numeric(d$INDELVAR)); st<-suppressWarnings(as.numeric(d[[tool]]))
    k<-!is.na(si)&!is.na(st)&!is.na(d$y); vid<-d$variant_id[k]; yk<-d$y[k]; a<-orient(si[k],yk); b<-orient(st[k],yk)
    ovm<-vid %in% ov; full<-d$variant_id %in% ov
    r1<-roc(yk[!ovm],a[!ovm],quiet=T,levels=c(0,1),direction="<")
    r2<-roc(yk[!ovm],b[!ovm],quiet=T,levels=c(0,1),direction="<")
    p<-roc.test(r1,r2,method="delong",paired=TRUE)$p.value
    data.table(Tool=tool,Class=ifelse(cls=="all","del+ins",cls),
      `Training overlap %`=sprintf(fmt,100*sum(full)/nrow(d)),
      `Leakage-free n`=sum(!ovm),`INDELVAR AUROC (LF)`=sprintf("%.3f",as.numeric(auc(r1))),
      `Tool AUROC (LF)`=sprintf("%.3f",as.numeric(auc(r2))),`Leakage-free DeLong p`=signif(p,2))}
  S5leak<-rbindlist(list(lfrow("SHINE","deletion"),lfrow("PON-Del","deletion"),
                         lfrow("MutPred-Indel","deletion"),
                         lfrow("IndeLLM","all"),lfrow("FATHMM-indel","all","%.2f")))
  # footnote stats derived from the computed leakage block.
  fa_pct <- S5leak[Tool=="FATHMM-indel"]$`Training overlap %`
  mp_ov  <- readLines(file.path(LK,"mutpredindel_train_overlap_ids.txt")); mp_ov<-mp_ov[nzchar(mp_ov)]
  mp_pd  <- m[indel_type=="deletion" & y==1]
  mp_nov <- sum(mp_pd$variant_id %in% mp_ov); mp_ntot <- nrow(mp_pd)
  sh_ov  <- readLines(file.path(LK,"shine_train_overlap_ids.txt")); sh_ov<-sh_ov[nzchar(sh_ov)]
  sh_bd  <- m[indel_type=="deletion" & y==0]
  sh_bpct<- sprintf("%.1f",100*sum(sh_bd$variant_id %in% sh_ov)/nrow(sh_bd))
  sh_pct <- S5leak[Tool=="SHINE"]$`Training overlap %`
  foot <- paste0("* training-set overlap with the test set; FATHMM-indel overlap (",fa_pct,"%) is negligible ",
             "and left unmarked. SHINE, PON-Del, and MutPred-Indel are evaluated on deletions (PON-Del is ",
             "deletion-only; SHINE and MutPred-Indel lacked adequate insertion coverage in this test set); ",
             "overlap % is over all test variants of that class. MutPred-Indel is trained on HGMD (pathogenic) ",
             "and rare gnomAD (neutral); its overlap is entirely on the pathogenic side (",round(100*mp_nov/mp_ntot),
             "%, ",format(mp_nov,big.mark=","),"/",format(mp_ntot,big.mark=",")," test ",
             "pathogenic deletions, by HGMD-accession match) with 0% on the common-benign side, so its full-",
             "coverage AUROC is the most leakage-inflated of any comparator. SHINE's overlap, by contrast, is ",
             "concentrated on the benign side (",sh_bpct,"% of test benign deletions vs ",sh_pct,"% of all test ",
             "deletions), inflating its full-coverage comparison; it is removed in the leakage-free reanalysis.")
  # Table S4: benchmark table only (asterisks explained in Table S5).
  s4<-file.path(TAB,"table_S4_test_set_benchmark.csv"); fwrite(S5bench,s4)
  cat("\n* training-set overlap with the test set; see Table S5 (leakage sensitivity).\n",file=s4,append=TRUE)
  write.xlsx(S5bench,file.path(TAB,"table_S4_test_set_benchmark.xlsx"))
  # Table S5: leakage / overlap sensitivity + full footnote.
  s5<-file.path(TAB,"table_S5_leakage_sensitivity.csv"); fwrite(S5leak,s5)
  cat(paste0("\nLeakage sensitivity: INDELVAR vs asterisk-marked tools after removing each tool's training overlap.\n",
             foot,"\n"),file=s5,append=TRUE)
  write.xlsx(S5leak,file.path(TAB,"table_S5_leakage_sensitivity.xlsx"))
  say("  wrote table_S4_test_set_benchmark + table_S5_leakage_sensitivity (.csv/.xlsx; %d benchmark, %d leakage rows)",nrow(S5bench),nrow(S5leak))
}

## ---------- Table S7: gnomAD in-frame-indel distribution across evidence bins ----------
sc  <- fread(file.path(PAPER,"gnomad_distribution_scored.tsv"))
ord <- c("PP3_+4","PP3_+3","PP3_+2","PP3_+1","No-evidence","BP4_-1","BP4_-2","BP4_-3","BP4_-4")
lbl <- function(p) fifelse(is.na(p)|p==0,"No-evidence",fifelse(p>0,paste0("PP3_+",p),paste0("BP4_",p)))
sc[, bin := lbl(acmg_points)]
S7 <- rbindlist(lapply(c("deletion","insertion"), function(ty){
  m <- sc[indel_type==ty]; tb <- table(factor(m$bin,levels=ord)); pc <- round(100*tb/sum(tb),2)
  say("  %s: gnomAD n=%d | PP3 %.1f%% (Strong+4 %.2f%%) | BP4 %.1f%%",
      ty, sum(tb), sum(pc[grep("PP3",ord)]), pc["PP3_+4"], sum(pc[grep("BP4",ord)]))
  data.table(type=ty,bin=ord,n=as.integer(tb),pct=as.numeric(pc)) }))
wr(S7,"table_S7_gnomad_distribution")
