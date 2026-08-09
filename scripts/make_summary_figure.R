args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop('Usage: Rscript make_summary_figure.R <results_dir> <output_prefix>')
resdir <- args[1]; prefix <- args[2]
m <- read.csv(file.path(resdir,'gse163256_reanalysis_models.csv'))
s <- read.csv(file.path(resdir,'signature_sensitivity_models.csv'))
pick <- function(df, outcome, model=NULL, signature=NULL, predictor=NULL) {
  z <- df[df$outcome==outcome,,drop=FALSE]
  if (!is.null(model)) z <- z[z$model==model,,drop=FALSE]
  if (!is.null(signature)) z <- z[z$signature==signature,,drop=FALSE]
  if (!is.null(predictor)) z <- z[z$predictor==predictor,,drop=FALSE]
  z[1,]
}
inf <- pick(m,'Inflamm','overall_mixed'); hyp <- pick(m,'Hypoxia','overall_mixed')
infp <- pick(s,'Inflamm',signature='all_overlap_pruned',predictor='Gly-minus-Ox')
hypp <- pick(s,'Hypoxia',signature='all_overlap_pruned',predictor='Gly-minus-Ox')
cols <- c(dark='#12344D',blue='#2C7FB8',teal='#2A9D8F',orange='#E69F00',red='#CC4C4C',purple='#7A5195',grey='#6B778D')
draw <- function() {
  par(mar=c(0,0,0,0),xaxs='i',yaxs='i'); plot.new(); plot.window(c(0,1),c(0,1))
  boxr <- function(x,y,w,h,fill,border){rect(x-w/2,y-h/2,x+w/2,y+h/2,col=fill,border=border,lwd=1.5)}
  txt <- function(x,y,z,cex=0.9,font=1,col='black',adj=.5) text(x,y,z,cex=cex,font=font,col=col,adj=adj)
  txt(.04,.95,'Evidence-aligned immunometabolic analysis in human intracerebral hemorrhage',1.55,2,cols['dark'],0);
  txt(.04,.91,'Verified data types, corrected longitudinal unit of analysis, sensitivity tests, and interpretation boundary',.9,1,cols['grey'],0)
  txt(.08,.84,'a   Verified datasets',1.15,2,cols['dark']); txt(.42,.84,'b   Corrected primary models',1.15,2,cols['dark']); txt(.78,.84,'c   Interpretation boundary',1.15,2,cols['dark'])
  boxr(.20,.69,.30,.19,'#EAF3F8',cols['blue']); txt(.20,.75,'GSE163256',1.2,2,cols['blue']); txt(.20,.68,'Longitudinal sorted-cell transcriptomics\nCD14+ monocytes/macrophages and neutrophils\n21 ICH patients + 5 healthy controls',.78)
  boxr(.20,.47,.30,.16,'#EAF7F4',cols['teal']); txt(.20,.52,'Locked monocyte analysis subset',1.0,2,cols['teal']); txt(.20,.45,'136 unique ICH observations\n20 patients; days 1-6\ntechnical duplicates collapsed',.82)
  boxr(.20,.24,.30,.16,'#FFF7E6',cols['orange']); txt(.20,.29,'GSE166638',1.0,2,cols['orange']); txt(.20,.22,'scRNA-seq from one longitudinally sampled patient\nexploratory cellular triangulation only',.76)
  boxr(.56,.68,.34,.20,'#F1F4F6',cols['dark']); txt(.56,.76,'Patient-level mixed models',1.0,2,cols['dark']); txt(.56,.68,sprintf('Inflammation: beta = %.3f (95%% CI %.3f to %.3f)\nHypoxia: beta = %.3f (95%% CI %.3f to %.3f)',inf$beta,inf$ci_low,inf$ci_high,hyp$beta,hyp$ci_low,hyp$ci_high),.79); txt(.56,.60,'outcome ~ contrast + compartment + day + (1 | patient)',.66,1,cols['grey'])
  boxr(.56,.43,.34,.18,'#F7F3FA',cols['purple']); txt(.56,.49,'Overlap-pruned sensitivity',1.0,2,cols['purple']); txt(.56,.42,sprintf('Inflammation: beta = %.3f (95%% CI %.3f to %.3f)\nHypoxia: beta = %.3f (95%% CI %.3f to %.3f)',infp$beta,infp$ci_low,infp$ci_high,hypp$beta,hypp$ci_low,hypp$ci_high),.74); txt(.56,.36,'All cross-signature overlaps removed',.63,1,cols['grey'])
  boxr(.87,.66,.21,.25,'#FFF1F1',cols['red']); txt(.87,.75,'Supported',1.05,2,cols['red']); txt(.87,.64,'Expression-derived metabolic scores\nco-vary with hypoxia and inflammation\nin longitudinal sorted monocyte profiles',.70)
  boxr(.87,.35,.21,.23,'#F5F5F5',cols['grey']); txt(.87,.43,'Not established',1.05,2,cols['grey']); txt(.87,.33,'Metabolic flux or causal niche effects\npopulation-level scRNA-seq replication\nor therapeutic actionability',.70)
  segments(.04,.10,.96,.10,col='#D0D7DE'); txt(.50,.065,'CellChat and virtual-knockout results from GSE166638 are hypothesis-generating only.',.78,2,cols['dark'])
}
dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE)
pdf(paste0(prefix,'.pdf'),width=11,height=6.8,useDingbats=FALSE); draw(); dev.off()
svg(paste0(prefix,'.svg'),width=11,height=6.8); draw(); dev.off()
png(paste0(prefix,'.png'),width=3300,height=2040,res=300); draw(); dev.off()
tiff(paste0(prefix,'.tiff'),width=3300,height=2040,res=300,compression='lzw'); draw(); dev.off()
