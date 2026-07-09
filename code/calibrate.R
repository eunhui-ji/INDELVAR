#!/usr/bin/env Rscript
## calibrate.R: ACMG evidence-strength cutoffs from the raw OOF predictions.
## (A) per-type ACMG/AMP evidence cutoffs -> evidence_cutoffs.tsv (Table 1)
## (B) local-LR+(s) 1-sided-95% cutoffs -> cutoffs_local_1sided95.tsv
## Usage: Rscript calibrate.R
suppressPackageStartupMessages({library(data.table)})
S   <- Sys.getenv("INDELVAR_ROOT",".")
M   <- file.path(S, "model")                                # read inputs
OUT <- Sys.getenv("INDELVAR_CUTOFFS_OUT", M)                 # write cutoffs (default = M)
say <- function(...) cat(sprintf(...),"\n")
need <- function(p) { if (!file.exists(p)) stop(sprintf("missing input: %s", p), call. = FALSE); p }

## (A) per-type evidence cutoffs: negative=ClinVar-B (OOF), gnomAD=window(>=3%)+prior;
## priors del .046/ins .008.
set.seed(1)
oof<-fread(need(file.path(S,"dataset","train_data","train_oof_predictions.tsv"))); setnames(oof,grep("indelvar|score|prob",names(oof),value=T)[1],"sc")
tr<-fread(need(file.path(S,"dataset","train_data","train_set_coords.tsv")))[,.(variant_id,indel_type)]; oof<-merge(oof,tr,by="variant_id")
abd<-fread(file.path(S,"dataset", "calibration", "gnomad_benign_scored.tsv")); abd[,it:=ifelse(indel_type=="deletion","deletion","insertion")]
source(file.path(S,"code","calib_params.R")); PARAM<-CALIB_PARAM   # alpha/OP
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
  allsc<-c(x,g); thrs<-sort(unique(c(allsc,floor(min(allsc)),ceiling(max(allsc)))),decreasing=TRUE); thrs_b<-rev(thrs)  # raw-unique thresholds
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

## (B) local-LR+(s) 1-sided-95% cutoffs. Targets = per-type ACMG target LRs
## (pooled = generic evidence-strength). iso_cutoff retained as NA for schema only.
set.seed(42)
E<-0.05; NBOOT<-10000; PC<-0.5
TARGETS_BY<-list(
  deletion =c(Supporting=2.39, Moderate=5.69,  Moderate3=13.59, Strong=32.42),
  insertion=c(Supporting=3.25, Moderate=10.55, Moderate3=34.28, Strong=111.34),
  pooled   =c(Supporting=2.08, Moderate=4.33,  Moderate3=9.00,  Strong=18.7))
GRID<-seq(0,1,by=0.005)

oof <- fread(need(file.path(S,"dataset","train_data","train_oof_predictions.tsv")))   # indelvar (raw), y, variant_id
tr2 <- fread(need(file.path(S,"dataset","train_data","train_set_coords.tsv")))[,.(variant_id, indel_type)]
oof <- merge(oof, tr2, by="variant_id", all.x=TRUE)
oof <- oof[!is.na(indelvar) & !is.na(y)]

local_lr<-function(x,ispath,idx){
  xb<-x[idx];pb<-ispath[idx];n1<-sum(pb);n0<-length(pb)-n1
  o<-order(xb);xs<-xb[o];ps<-pb[o];cumP<-cumsum(ps);cumN<-cumsum(1-ps)
  lo<-GRID-E;hi<-GRID+E
  upP<-cumP[pmax(findInterval(hi,xs),1)];upP[findInterval(hi,xs)<1]<-0
  dnP<-cumP[pmax(findInterval(lo,xs),1)];dnP[findInterval(lo,xs)<1]<-0
  upN<-cumN[pmax(findInterval(hi,xs),1)];upN[findInterval(hi,xs)<1]<-0
  dnN<-cumN[pmax(findInterval(lo,xs),1)];dnN[findInterval(lo,xs)<1]<-0
  c1<-upP-dnP;c0<-upN-dnN
  ((c1+PC)/(n1+PC))/((c0+PC)/(n0+PC))
}
derive<-function(sub){
  d<-if(sub=="pooled")oof else oof[indel_type==sub]
  TARGETS<-TARGETS_BY[[sub]]
  x<-d$indelvar;ispath<-d$y;n<-length(x)
  BM<-matrix(NA_real_,NBOOT,length(GRID))
  for(b in 1:NBOOT) BM[b,]<-local_lr(x,ispath,sample.int(n,n,replace=TRUE))
  LB<-apply(BM,2,quantile,0.05,na.rm=TRUE);UB<-apply(BM,2,quantile,0.95,na.rm=TRUE)
  res<-list()
  for(nm in names(TARGETS)){tg<-TARGETS[[nm]]
    okP<-LB>=tg&is.finite(LB);tP<-NA_real_;run<-TRUE
    for(j in length(GRID):1){if(!okP[j])run<-FALSE;if(run)tP<-GRID[j]}
    okB<-UB<=1/tg&is.finite(UB);tB<-NA_real_;run<-TRUE
    for(j in 1:length(GRID)){if(!okB[j])run<-FALSE;if(run)tB<-GRID[j]}
    res[[paste0("PP3_",nm)]]<-tP;res[[paste0("BP4_",nm)]]<-tB}
  res
}
out<-list()
for(sub in c("deletion","insertion","pooled")){r<-derive(sub); TG<-TARGETS_BY[[sub]]
  for(ev in names(r)){raw<-r[[ev]]; nm<-sub("PP3_|BP4_","",ev)
    tgt<-if(grepl("^PP3",ev)) TG[[nm]] else 1/TG[[nm]]
    out[[length(out)+1]]<-data.table(subset=ifelse(sub=="deletion","del",ifelse(sub=="insertion","ins","pooled")),
      evidence=ev, target=round(tgt,4), raw_cutoff=raw, iso_cutoff=NA_real_, reachable=!is.na(raw))}
  cat(sprintf("[R2] %s reachable=%d/8\n",sub,sum(sapply(r,function(z)!is.na(z)))))}
DT<-rbindlist(out)
DT[,raw_cutoff:=round(raw_cutoff,4)][,iso_cutoff:=round(iso_cutoff,4)]
setcolorder(DT,c("subset","evidence","target","raw_cutoff","iso_cutoff","reachable"))
print(DT)
fwrite(DT, file.path(OUT,"cutoffs_local_1sided95.tsv"), sep="\t")
cat("[R2] saved cutoffs_local_1sided95.tsv\n")
