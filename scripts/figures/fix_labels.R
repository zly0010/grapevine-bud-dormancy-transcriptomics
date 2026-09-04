options(stringsAsFactors=FALSE)
# Package library paths are supplied by the user environment.
stopifnot(requireNamespace('ragg'),requireNamespace('jsonlite'),requireNamespace('systemfonts'))
library(grid)
root <- 'outputs/figure_labels'
src <- 'outputs/advanced/results/29_figure_layout_revision_2026-08-07/figures'
# Figure contract: retain every original data pixel, scale, panel, color and number.
# Quantitative layouts and scientific conclusions are frozen. Only four labels change.
# R-only raster decoding, label rendering, compositing, PNG export and pixel QA.
read_png <- function(path) {
  b <- readBin(path,'raw',n=file.info(path)$size)
  uint <- function(i) sum(as.integer(b[i+(0:3)])*256^(3:0))
  stopifnot(identical(as.integer(b[1:8]),c(137L,80L,78L,71L,13L,10L,26L,10L)))
  at <- 9; z <- raw()
  while(at < length(b)) {
    n <- uint(at); kind <- rawToChar(b[(at+4):(at+7)])
    dat <- b[(at+8):(at+7+n)]
    if(kind=='IHDR') {w<-uint(at+8);h<-uint(at+12); depth<-as.integer(dat[9]);type<-as.integer(dat[10]);stopifnot(depth==8, type %in% c(2,6),dat[13]==0)}
    if(kind=='IDAT') z<-c(z,dat)
    if(kind=='IEND') break
    at<-at+n+12
  }
  ch<-if(type==2) 3L else 4L; stride<-w*ch
  bytes<-as.integer(memDecompress(z,type='gzip')); stopifnot(length(bytes)==h*(stride+1))
  out<-matrix(0L,h,stride); prev<-integer(stride)
  for(y in seq_len(h)) {
    offset<-(y-1)*(stride+1); filter<-bytes[offset+1]; a<-bytes[offset+1+seq_len(stride)]
    if(filter==1) for(k in seq_len(ch)) {ii<-seq.int(k,stride,ch);a[ii]<-cumsum(a[ii])%%256}
    if(filter==2) a<-(a+prev)%%256
    if(filter %in% c(3,4)) for(x in seq_len(stride)) {
      left<-if(x>ch) a[x-ch] else 0L; up<-prev[x]; ul<-if(x>ch) prev[x-ch] else 0L
      if(filter==3) p<-floor((left+up)/2) else {p<-left+up-ul;pa<-abs(p-left);pb<-abs(p-up);pc<-abs(p-ul);p<-if(pa<=pb && pa<=pc) left else if(pb<=pc) up else ul}
      a[x]<-(a[x]+p)%%256
    }
    stopifnot(filter %in% 0:4);out[y,]<-a;prev<-a
  }
  if(ch==4) stopifnot(all(out[,seq(4,stride,4)]==255))
  colors<-rgb(out[,seq(1,stride,ch)],out[,seq(2,stride,ch)],out[,seq(3,stride,ch)],maxColorValue=255)
  matrix(colors,h,w)
}
bbox <- function(m) {p<-which(m!='#FFFFFF',arr.ind=TRUE);stopifnot(nrow(p)>0);c(min(p[,2]),min(p[,1]),max(p[,2]),max(p[,1]))}
render_text <- function(label,size,bold) {
  capture<-ragg::agg_capture(width=1600,height=130,res=300,background='white')
  grid.text(label,x=unit(10,'native'),y=.5,just=c('left','centre'),gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold) 'bold' else 'plain',col='#202124'),vp=viewport(xscale=c(0,1600)))
  m<-capture();dev.off();palette<-unique(as.vector(m));rgbvals<-col2rgb(palette);hex<-rgb(rgbvals[1,],rgbvals[2,],rgbvals[3,],maxColorValue=255);m<-matrix(hex[match(m,palette)],nrow(m),ncol(m));bb<-bbox(m);m[bb[2]:bb[4],bb[1]:bb[3],drop=FALSE]
}
specs <- list(
  list(id=3,stem='Figure_3_module_stability_and_trajectories_R',sw=510.24,sh=425.20,
       labels=list(list(old='Leave-one-variety module recovery',new='Leave-one-genotype module recovery',box=c(294,7,449,18.8),size=8,bold=TRUE,align='left'),list(old='Omitted variety',new='Omitted genotype',box=c(338,165.5,412,179),size=7.5,bold=FALSE,align='center'))),
  list(id=6,stem='figure_6_virtual_perturbation',sw=510.24,sh=496.06,
       labels=list(list(old='Four leave-one-variety refits',new='Four leave-one-genotype refits',box=c(204,13,296.5,24),size=6.7,bold=FALSE,align='left'),list(old='Held-out variety',new='Held-out genotype',box=c(374,476.5,447,488),size=7.5,bold=FALSE,align='center')))
)
logs<-list()
for(s in specs) {
  input<-file.path(src,paste0(s$stem,'.png')); original<-read_png(input); updated<-original
  cat('Decoded figure',s$id,dim(original),'\n');flush.console()
  height<-nrow(original);width<-ncol(original);allowed<-matrix(FALSE,height,width)
  svg<-readLines(file.path(src,paste0(s$stem,'.svg')),warn=FALSE)
  # Use actual frozen viewBox dimensions, not inferred physical sizing.
  head<-svg[grepl('<svg ',svg)][1]; vb<-strsplit(sub(".*viewBox='([^']+)'.*",'\\1',head),' ')[[1]]
  s$sw<-as.numeric(vb[3]);s$sh<-as.numeric(vb[4])
  entries<-list()
  for(l in s$labels) {
    b<-c(floor(l$box[1]/s$sw*width),floor(l$box[2]/s$sh*height),ceiling(l$box[3]/s$sw*width),ceiling(l$box[4]/s$sh*height))
    region<-original[b[2]:b[4],b[1]:b[3],drop=FALSE]; oldbb<-bbox(region)
    oldpatch<-render_text(l$old,l$size,l$bold); patch<-render_text(l$new,l$size,l$bold)
    # Check font size/weight against the original measured glyph bounds.
    cat(l$old,'original',oldbb[3]-oldbb[1]+1,oldbb[4]-oldbb[2]+1,'rendered',ncol(oldpatch),nrow(oldpatch),'\n')
    stopifnot(abs(nrow(oldpatch)-(oldbb[4]-oldbb[2]+1))<=2,abs(ncol(oldpatch)-(oldbb[3]-oldbb[1]+1))<=4)
    x<-if(l$align=='left') b[1]+oldbb[1]-1 else round(b[1]+(oldbb[1]+oldbb[3])/2-1-(ncol(patch)-1)/2)
    y<-b[2]+oldbb[2]-1
    stopifnot(x>=b[1],x+ncol(patch)-1<=b[3],y+nrow(patch)-1<=b[4])
    updated[b[2]:b[4],b[1]:b[3]]<-'#FFFFFF'
    updated[y+seq_len(nrow(patch))-1,x+seq_len(ncol(patch))-1]<-patch
    allowed[b[2]:b[4],b[1]:b[3]]<-TRUE
    idx<-grep(l$old,svg,fixed=TRUE);stopifnot(length(idx)==1)
    ratio<-systemfonts::string_width(l$new,family='Arial',bold=l$bold,size=l$size)/systemfonts::string_width(l$old,family='Arial',bold=l$bold,size=l$size)
    oldlen<-as.numeric(sub(".*textLength='([0-9.]+)px'.*",'\\1',svg[idx]))
    svg[idx]<-sub(l$old,l$new,svg[idx],fixed=TRUE)
    svg[idx]<-sub("textLength='[0-9.]+px'",sprintf("textLength='%.2fpx'",oldlen*ratio),svg[idx])
    entries[[length(entries)+1]]<-list(old=l$old,new=l$new,allowed_pixel_rectangle_1based=b,old_glyph_dimensions=c(oldbb[3]-oldbb[1]+1,oldbb[4]-oldbb[2]+1),rendered_old_dimensions=c(ncol(oldpatch),nrow(oldpatch)))
  }
  stopifnot(all(updated[!allowed]==original[!allowed]),!any(grepl('variety',svg,ignore.case=TRUE)))
  output<-file.path(root,paste0('figure',s$id,'_label_fixed.png'))
  ragg::agg_png(output,width=width,height=height,res=300,background='white')
  grid.raster(as.raster(updated),width=1,height=1,interpolate=FALSE);dev.off()
  verified<-read_png(output);cat('Export pixel differences:',sum(updated!=verified),'\n'); print(head(cbind(updated[updated!=verified],verified[updated!=verified])));stopifnot(all(updated==verified))
  writeLines(svg,file.path(root,paste0('figure',s$id,'_label_fixed.svg')),useBytes=TRUE)
  logs[[length(logs)+1]]<-list(figure=s$id,original=input,updated=output,width=width,height=height,labels=entries,changed_pixels=sum(updated!=original),changed_pixels_outside_label_regions=sum(updated[!allowed]!=original[!allowed]),data_geometry_scales_colors_numbers_unchanged=TRUE)
  cat('Verified figure',s$id,'outside-label differences = 0\n');flush.console()
}
jsonlite::write_json(logs,file.path(root,'qa/FIGURE_LABEL_VERIFICATION.json'),pretty=TRUE,auto_unbox=TRUE)
