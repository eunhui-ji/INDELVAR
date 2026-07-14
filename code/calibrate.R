#!/usr/bin/env Rscript
## calibrate.R: ACMG evidence-strength cutoffs and LR targets from the raw OOF predictions.
suppressPackageStartupMessages({library(data.table)})
S   <- Sys.getenv("INDELVAR_ROOT",".")
M   <- file.path(S, "model")
OUT <- Sys.getenv("INDELVAR_CUTOFFS_OUT", M)
say <- function(...) cat(sprintf(...),"\n")
need <- function(p) { if (!file.exists(p)) stop(sprintf("missing input: %s", p), call. = FALSE); p }

set.seed(1)
oof<-fread(need(file.path(S,"dataset","train_data","train_oof_predictions.tsv"))); setnames(oof,grep("indelvar|score|prob",names(oof),value=T)[1],"sc")
tr<-fread(need(file.path(S,"dataset","train_data","train_set_coords.tsv")))[,.(variant_id,indel_type)]; oof<-merge(oof,tr,by="variant_id")
abd<-fread(file.path(S,"dataset", "calibration", "gnomad_calibration_scored.tsv")); abd[,it:=ifelse(indel_type=="deletion","deletion","insertion")]
source(file.path(S,"code","calib_params.R")); PARAM<-CALIB_PARAM
minpoints<-100; gft<-0.03; increment<-0.001; B<-10000; disc<-0.05; HGRID<-seq(0,1,by=increment)
EV<-c("PP3_Supporting","PP3_Moderate","PP3_Moderate3","PP3_Strong","BP4_Supporting","BP4_Moderate","BP4_Moderate3","BP4_Strong")
PTS<-c(1,2,3,4,-1,-2,-3,-4)
local_post<-function(x,y,g,thrs,w){ maxt<-max(thrs);mint<-min(thrs);Ng<-length(g)
  xp<-sort(x[y==1]); xn<-sort(x[y==0]); gs<-sort(g)
  cnt<-function(v,lo,hi) findInterval(hi,v)-findInterval(lo-1e-12,v)
  sapply(thrs,function(t){ ok<-function(h){lo<-t-h;hi<-t+h
      cc<-if(hi>maxt)(maxt-lo)/(hi-lo) else if(lo<mint)(hi-mint)/(hi-lo) else 1; if(cc<=0)cc<-1e-9
      (cnt(xp,lo,hi)+cnt(xn,lo,hi) >= cc*minpoints) && (cnt(gs,lo,hi) >= gft*cc*Ng)}
    lo<-1L;hi<-length(HGRID); if(!ok(HGRID[hi])) return(NA_real_)
    while(lo<hi){m<-(lo+hi)%/%2L; if(ok(HGRID[m])) hi<-m else lo<-m+1L}
    h<-HGRID[lo]; pos<-cnt(xp,t-h,t+h); neg<-cnt(xn,t-h,t+h); pos/(pos+w*neg)}) }
get_thr<-function(post,thrs,Post) sapply(seq_along(Post),function(j){
  ind<-suppressWarnings(min(which(post<Post[j]))-1); if(is.finite(ind)&&ind>0) thrs[ind] else NA_real_})
discf<-function(Mb,Post,B,d,type){ sapply(1:ncol(Mb),function(j){inv<-sum(is.na(Mb[,j]))
  if(inv> d*B) NA_real_ else {t<-sort(Mb[!is.na(Mb[,j]),j],decreasing=(type=="path")); t[floor(d*B)-inv+1]}})}
rows<-list()
for(tt in c("deletion","insertion")){ pp<-PARAM[[tt]]; alpha<-pp$alpha
  d<-oof[indel_type==tt & !is.na(sc) & !is.na(y)]; x<-d$sc; yv<-d$y
  g<-abd[it==tt & !is.na(indelvar_score),indelvar_score]
  w<-(1-alpha)*sum(yv==1)/(sum(yv==0)*alpha)
  LR<-pp$OP^((1:4)/8); Post_p<-LR*alpha/(LR*alpha+(1-alpha)); Post_b<-LR*(1-alpha)/((LR-1)*(1-alpha)+1)
  allsc<-c(x,g); thrs<-sort(unique(c(allsc,floor(min(allsc)),ceiling(max(allsc)))),decreasing=TRUE); thrs_b<-rev(thrs)
  BMp<-matrix(NA_real_,B,4); BMb<-matrix(NA_real_,B,4)
  for(b in 1:B){qx<-sample.int(length(x),replace=TRUE);qg<-sample.int(length(g),replace=TRUE)
    post<-local_post(x[qx],yv[qx],g[qg],thrs,w); postb<-1-rev(post)
    BMp[b,]<-get_thr(post,thrs,Post_p); BMb[b,]<-get_thr(postb,thrs_b,Post_b)}
  DTp<-discf(BMp,Post_p,B,disc,"path"); DTb<-discf(BMb,Post_b,B,disc,"benign")
  g3<-if(tt=="deletion")"del" else "ins"
  say("%s w=%.1f | PP3 disc=%s | BP4 disc=%s",tt,w,paste(round(DTp,3),collapse="/"),paste(round(DTb,3),collapse="/"))
  for(k in 1:4) rows[[length(rows)+1]]<-data.table(subset=g3,evidence=EV[k],points=PTS[k],cutoff=DTp[k],reached=!is.na(DTp[k]))
  for(k in 1:4) rows[[length(rows)+1]]<-data.table(subset=g3,evidence=EV[k+4],points=PTS[k+4],cutoff=DTb[k],reached=!is.na(DTb[k]))
}
res<-rbindlist(rows); fwrite(res,file.path(OUT,"evidence_cutoffs.tsv"),sep="\t")
say("\nwrote evidence_cutoffs.tsv"); print(res)

TARGETS_BY<-list(
  deletion =c(Supporting=2.39, Moderate=5.69,  Moderate3=13.59, Strong=32.42),
  insertion=c(Supporting=3.25, Moderate=10.55, Moderate3=34.28, Strong=111.34),
  pooled   =c(Supporting=2.08, Moderate=4.33,  Moderate3=9.00,  Strong=18.7))
trows<-list()
for(sub in names(TARGETS_BY)){ TG<-TARGETS_BY[[sub]]
  g3<-ifelse(sub=="deletion","del",ifelse(sub=="insertion","ins","pooled"))
  for(nm in names(TG)){
    trows[[length(trows)+1]]<-data.table(subset=g3,evidence=paste0("PP3_",nm),target=round(TG[[nm]],4))
    trows[[length(trows)+1]]<-data.table(subset=g3,evidence=paste0("BP4_",nm),target=round(1/TG[[nm]],4))}}
TGT<-rbindlist(trows)
print(TGT)
fwrite(TGT, file.path(OUT,"evidence_strength_targets.tsv"), sep="\t")
cat("saved evidence_strength_targets.tsv\n")
