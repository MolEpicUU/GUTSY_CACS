rm(list=ls())

#set the working director
set.seed(123)
input.path="./1_simulations/"
output.path="./Ext.Data.Fig/"
dir.create(output.path, showWarnings = FALSE)
ncores=1

#load libraries
library(rio)
library(ggplot2)
library(gridExtra)

#import results
res_spearman=import(paste(input.path,"res_fastlm.fun_CASC_log1pspearman_1batch.tsv",sep=""))
res_lm=import(paste(input.path,"res_fastlm.fun_CASC_log1p_1batch.tsv",sep=""))
res_lm.clr=import(paste(input.path,"res_fastlm.fun_CASC_clr_1batch.tsv",sep=""))
res_ordinal=import(paste(input.path,"res_ordinal.fun_CASC_log1p_1batch.tsv",sep=""))
res_two.sample=import(paste(input.path,"res_two.stage.log.lm.fun_log1p_1batch.tsv",sep=""))
res_nb.clr=import(paste(input.path,"nb_clr_1batch.tsv",sep=""))
res_nb=import(paste(input.path,"nb_log1p_1batch.tsv",sep=""))
res_hurdle.nb=import(paste(input.path,"res_hurdle.nb.fun_log1p_1batch.tsv",sep=""))
res_boot.lm=import(paste(input.path,"res_fastboot.fun_CASC_log1p_1batch.tsv",sep=""))
res_boot.lm.residual=import(paste(input.path,"res_fastboot_resi.fun_CASC_log1p_1batch.tsv",sep=""))
res_boot.nb=import(paste(input.path,"res_nb.boot.log1p_1batch.tsv",sep=""))
res_rob=import(paste(input.path,"res_fastrobust.fun_CASC_log1p_1batch.tsv",sep=""))

#load functions
source("/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/Demo/functions/estlambda.fun.R")


#remove na
res_spearman=na.omit(res_spearman)
res_lm=na.omit(res_lm)
res_lm.clr=na.omit(res_lm.clr)
res_ordinal=na.omit(res_ordinal)
res_two.sample=na.omit(res_two.sample)
res_nb=na.omit(res_nb)
res_nb.clr=na.omit(res_nb.clr)
res_hurdle.nb=na.omit(res_hurdle.nb)
res_boot.lm=na.omit(res_boot.lm)
res_boot.nb=na.omit(res_boot.nb)
res_boot.lm.residual=na.omit(res_boot.lm.residual)
res_rob=na.omit(res_rob)


for(i in grep("res",ls(),value=T))
{
  try(eval(parse(text=paste(i,"$p.value[",i,"$p.value<1e-200]=1e-200",sep=""))))
}


#order the pvalues and estimate the expected pvalues for the different methods
a=sort(res_spearman$p.value)
b=sort(res_lm$p.value)
ae=sort(ppoints(res_spearman$p.value))
be=sort(ppoints(res_lm$p.value))
b.clr=sort(res_lm.clr$p.value)
be.clr=sort(ppoints(res_lm.clr$p.value))

c=sort(res_ordinal$p.value)
ce=sort(ppoints(res_ordinal$p.value))
d=sort(res_two.sample$p.value.zero)
de=sort(ppoints(res_two.sample$p.value.zero))
e=sort(res_two.sample$p.value.count)
ee=sort(ppoints(res_two.sample$p.value.count))
f=sort(res_nb$p.value)
ff=sort(ppoints(res_nb$p.value))
f.clr=sort(res_nb.clr$p.value)
ff.clr=sort(ppoints(res_nb.clr$p.value))
g=sort(res_hurdle.nb$p.value.zero)
gg=sort(ppoints(res_hurdle.nb$p.value.zero))
h=sort(res_hurdle.nb$p.value.count)
hh=sort(ppoints(res_hurdle.nb$p.value.count))
i=sort(res_boot.lm$p.value)
ii=sort(ppoints(res_boot.lm$p.value))
j=sort(res_boot.nb$p.value)
jj=sort(ppoints(res_boot.nb$p.value))
k=sort(res_boot.lm.residual$p.value)
kk=sort(ppoints(res_boot.lm.residual$p.value))
l=sort(res_rob$p.value)
ll=sort(ppoints(res_rob$p.value))

#create a dataframe with the ordered p.values and the expected p.values
a0=data.frame(a,ae)
b0=data.frame(b,be)
b0.clr=data.frame(b.clr,be.clr)
c0=data.frame(c,ce)
d0=data.frame(d,de)
e0=data.frame(e,ee)
f0=data.frame(f,ff)
f0.clr=data.frame(f.clr,ff.clr)
g0=data.frame(g,gg)
h0=data.frame(h,hh)
i0=data.frame(i,ii)
j0=data.frame(j,jj)
k0=data.frame(k,kk)
l0=data.frame(l,ll)

#identify the maximum number of p.values detected among all the methods 
max.c=max(nrow(a0),nrow(b0),nrow(b0.clr),
          nrow(c0),
          nrow(d0),
          nrow(e0),
          nrow(f0),
          nrow(f0.clr),
          nrow(g0),nrow(h0),
          nrow(i0),
          nrow(j0),
          nrow(k0),nrow(l0)
          )

#make all the dataframe the same length
a0=a0[1:max.c,]
b0=b0[1:max.c,]
b0.clr=b0.clr[1:max.c,]
c0=c0[1:max.c,]
d0=d0[1:max.c,]
e0=e0[1:max.c,]
f0=f0[1:max.c,]
f0.clr=f0.clr[1:max.c,]
g0=g0[1:max.c,]
h0=h0[1:max.c,]
i0=i0[1:max.c,]
j0=j0[1:max.c,]
k0=k0[1:max.c,]
l0=l0[1:max.c,]

#estimate the inflation lambda for each method
a.t1=round(estlambda(res_spearman$p.value)$estimate,2)
b.t1=round(estlambda(res_lm$p.value)$estimate,2)
b.t1.clr=round(estlambda(res_lm.clr$p.value)$estimate,2)
c.t1=round(estlambda(res_ordinal$p.value)$estimate,2)
d.t1=round(estlambda(res_two.sample$p.value.zero)$estimate,2)
e.t1=round(estlambda(res_two.sample$p.value.count)$estimate,2)
f.t1=round(estlambda(res_nb$p.value)$estimate,2)
f.t1.clr=round(estlambda(res_nb.clr$p.value)$estimate,2)
g.t1=round(estlambda(res_hurdle.nb$p.value.zero)$estimate,2)
h.t1=round(estlambda(res_hurdle.nb$p.value.count)$estimate,2)
i.t1=round(estlambda(res_boot.lm$p.value)$estimate,2)
j.t1=round(estlambda(res_boot.nb$p.value)$estimate,2)
k.t1=round(estlambda(res_boot.lm.residual$p.value)$estimate,2)
l.t1=round(estlambda(res_rob$p.value)$estimate,2)

#merge all the method p.values in one data.frame
aa=data.frame(a0,b0,b0.clr,
              c0,
              d0,e0,
              f0,
              f0.clr,g0,h0,
              i0,j0,
              k0,l0
              )

#create the qq-plot
g1=ggplot(data=aa,)+
  geom_abline()+
  geom_point(aes(x=-log10(ae),y=-log10(a),color="Spearman"),alpha=0.75)+
  geom_point(aes(x=-log10(be.clr),y=-log10(b.clr),color="lm.clr"),alpha=0.75)+
  geom_point(aes(x=-log10(de),y=-log10(d),color="two.sample.zero"),alpha=0.75)+
  geom_point(aes(x=-log10(ee),y=-log10(e),color="two.sample.count"),alpha=0.75)+
  geom_point(aes(x=-log10(ff),y=-log10(f),color="nb"),alpha=0.75)+
  geom_point(aes(x=-log10(ff.clr),y=-log10(f.clr),color="nb.clr"),alpha=0.75)+
  geom_point(aes(x=-log10(gg),y=-log10(g),color="hurdle.nb.zero"),alpha=0.75)+
  geom_point(aes(x=-log10(hh),y=-log10(h),color="hurdle.nb.count"),alpha=0.75)+
  geom_point(aes(x=-log10(ii),y=-log10(i),color="boot.lm"),alpha=0.75)+
  geom_point(aes(x=-log10(jj),y=-log10(j),color="boot_nb"),alpha=0.75)+
  geom_point(aes(x=-log10(ll),y=-log10(l),color="rob"),alpha=0.75)+
  geom_point(aes(x=-log10(ce),y=-log10(c),color="ordinal"),alpha=0.75)+
  geom_point(aes(x=-log10(be),y=-log10(b),color="lm.log1"),alpha=0.75)+
  geom_point(aes(x=-log10(kk),y=-log10(k),color="boot.lm.resi"),alpha=0.75)+
  
  
  # geom_point(aes(x=-log10(mm),y=-log10(m),color="boot_lm_resi2"),alpha=0.75)+
  scale_colour_manual(name="models",
                      values = c("tan","turquoise2","yellow",
                                 "steelblue3",
                                 "lightpink","green","blue",
                                 "darkmagenta", "darkseagreen4","black",
                                 "grey","orangered",
                                  "#E69F00","#999999"
                                 ),
                      labels = c(
                        paste("boot.nb (lambda=",j.t1,")",sep=""),
                        paste("boot.lm (lambda=",i.t1,")",sep=""),
                        paste("boot.lm.resi (lambda=",k.t1,")",sep=""),
                        # paste("boot_lm_resi2 (lambda=",m.t1,")",sep=""),
                        paste("hurdle.nb.count (lambda=",h.t1,")",sep=""),
                        paste("hurdle.nb.zero (lambda=",g.t1,")",sep=""),
                        paste("lm.clr (lambda=",b.t1.clr,")",sep=""),
                        paste("lm.log1p (lambda=",b.t1,")",sep=""),
                        paste("nb (lambda=",f.t1,")",sep=""),
                        paste("nb.clr (lambda=",f.t1.clr,")",sep=""),
                        paste("ordinal (lambda=",c.t1,")",sep=""),
                        paste("robust.lm (lambda=",l.t1,")",sep=""),
                        paste("spearman (lambda=",a.t1,")",sep=""),
                        paste("two.stage.count (lambda=",e.t1,")",sep=""),
                        paste("two.stage.zero (lambda=",d.t1,")",sep="")
                        
                      ))+
  xlab("-log10(exp.p.value)")+
  ylab("-log10(obs.p.value)")+
  ylim(0,10)+
  # use a nice theme:
  theme_bw()+
  # here you tweak all the micro settings: 
  theme(panel.grid.major.x = element_line(color = 'gray80', linetype = 3), 
        panel.grid.major.y = element_line(color = 'gray80', linetype = 3), 
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size=20),  
        axis.title.x = element_text(size=20),  
        axis.text.y = element_text(size=20), #element_blank(),# this makes the horizontal Strep labels disappear. If you want to allow them back in remove the element_blank() and add this element_text(size=9)
        axis.text.x = element_text(size=20),
        axis.ticks.length=unit(0,"cm"),
        
        legend.position="right",
        legend.key.size = unit(0.5,"cm"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=17),
        
        strip.placement = "inside",  # outside
        strip.text =  element_text(size=20, face = 'bold'),  # replace this with element_blank() if you turn axis.text.y into element_text
        
        panel.border = element_rect(color = "gray95", fill = NA, size = 1),
        panel.spacing = unit(0, "lines"), # space between facets
        
        strip.background =element_blank(), 
        
        panel.grid.minor.y = element_blank()) 


#save the qqplot
ggsave(g1,file=paste(output.path,"qqplot_models_casctot.pdf",sep=""),width=10,height=10)
ggsave(g1,file=paste(output.path,"qqplot_models_casctot.svg",sep=""),width=10,height=10)

