## eccDNA鉴定
Circle-Map是一款用于全基因组测序，尤其是环状DNA富集测序（Circle-Seq）中的eccDNA的识别。识别流程如下：

[windows下载circle-map并使用](https://www.yuque.com/white-f5mnu/qhi4bs/nou4va9ebw6qhsoh)

## 统计eccDNA的长度并可视化
示例数据：[示例数据.zip](https://www.yuque.com/attachments/yuque/0/2023/zip/38459666/1698223105361-ad3618eb-82f1-412b-a924-c6c6a09d6a35.zip)

示例数据说明：两个bed文件。bed文件使用tab分隔符文本文件，可以使用excel和记事本打开

数据分析工具：R

### 代码
```bash
# 设置工作路径
setwd("D:\\HepG2\\代码\\srtp")

# 导入两个文件定义为df1和df2（文件需要在工作路径内）
df1 <- read.csv(file = "example1.bed",header = FALSE,sep = "\t")
df2 <- read.csv(file = "example2.bed",header = FALSE,sep = "\t")

# 计算eccDNA的长度
df1$Length <- df1$V3-df1$V2
df2$Length <- df2$V3-df2$V2

# 给每一个数据框一个分组，定为control列和dox列
df1$group <- "Control"
df2$group <- "DOX"

# 查看当前数据框的信息
head(df1)

# 按行合并数据框
df3 <- rbind(df1,df2)

# 安装和载入包
install.packages("ggplot2")
install.packages("ggprism")
library(ggplot2)
library(ggprism)

# 显示每个数据的行数，即eccDNA的个数
example1_total_num <- length(df1$Length)
example2_total_num <- length(df2$Length)
print(example1_total_num);print(example2_total_num)

# 计算小于3000bp的eccDNA占总数的比例
length(subset(df1,Length<3000)$Length)/example1_total_num
length(subset(df2,Length<3000)$Length)/example2_total_num
#分别是[1] 0.9739958;[1] 0.9765895

# 使用ggplot绘图
ggplot(data = df3,mapping = aes(x=Length,color=group))+
  # adjust调整图的平滑度，linewidth调整线的粗细
  geom_density(adjust=0.2,linewidth=1)+
  # 调整横坐标的范围
  scale_x_continuous(limits = c(0,3000))+
  # y轴和x轴的标签
  labs(x = "eccDNA size(bp)", y = "Density")+
  # 添加主题类型,参数依次为：调整字体的大小，线的粗细，字的粗细
  theme_prism(base_size = 14,base_line_size = 1,base_fontface = "plain")+
  # 线条的颜色
  scale_color_manual(values = c("blue","red"))

# 保存图片为pdf格式
# pdf是矢量图，放大不会模糊，在科研绘图中使用pdf保存是比较稳妥的
# 适量图和位图的区别：https://jingyan.baidu.com/article/54b6b9c0dbef682d583b4722.html
ggsave(filename = "length_distribution.pdf",device = "pdf",width = 8,height = 4)


# 累计频率图
ggplot(df3, aes(Length, colour = group)) +
  stat_ecdf(linewidth=1.2)+
  # y轴和x轴的标签
  labs(x = "eccDNA size(bp)", y = "Cumulative frequency(%)")+
  # 调整横坐标的范围
  scale_x_continuous(limits = c(0,3000),breaks=seq(0,3000,200))+
  scale_y_continuous(breaks=c(0,0.25,0.50,0.75,1),labels=c("0","25","50","75","100"))+
  # 更改主题
  theme_prism(base_size = 14,base_line_size = 1,
              base_fontface = "plain",axis_text_angle = 45)+
  # 线条的颜色
  scale_color_manual(values = c("blue","red"))
  
# 保存图片为pdf格式。
ggsave(filename = "Cumulative_frequency_of_length.pdf",device = "pdf",width = 8,height = 4)
```



### 结果
![eccDNA的长度分布](https://cdn.nlark.com/yuque/0/2023/png/38459666/1698226398383-b4b50af2-039f-4ea1-8745-a6e3e99168b8.png)

example1和example2的eccDNA的数量分别为：4307和12174。长度低于3000bp的百分比为97.40%和97.66%，说明主要是小eccDNA。eccDNA的分布呈现200bp左右的周期性。

![eccDNA的累积频率](https://cdn.nlark.com/yuque/0/2023/png/38459666/1698227841805-a8a50af7-d4b1-47dd-827d-9df066f67da1.png)

累计频率图显示：example2的eccDNA长度比example1的更长。

## eccDNA在染色体上的分布并可视化
### 方法1：karyoploteR包
1. 绘图代码

```bash
# 设置工作路径
setwd("D:\\HepG2\\代码\\srtp")

# 导入两个文件定义为df1和df2（文件需要在工作路径内）
df1 <- read.csv(file = "example1.bed",header = FALSE,sep = "\t")
df2 <- read.csv(file = "example2.bed",header = FALSE,sep = "\t")

# 如果BiocManager不存在就安装BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# 使用BiocManager安装"karyoploteR"和"GenomicRanges"包
BiocManager::install(c("karyoploteR","GenomicRanges"))

# 载入需要的包
library("karyoploteR")
library("GenomicRanges")

# 创建GRanges对象，后续kpPlotDensity需要这种格式的数据
gr_example1 <- GRanges(
  seqnames = df1$V1,
  ranges = IRanges(start = df1$V2,end = df1$V3),
  strand = "*"
)

gr_example2 <- GRanges(
  seqnames = df2$V1,
  ranges = IRanges(start = df2$V2,end = df2$V3),
  strand = "*"
)

# 显示gr_example1的头部
head(gr_example1)

# 画密度图
# 教程详见https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotDensity/PlotDensity.html
# 保存pdf()和dev.off()中间代码生成的图片，命名为Location of eccDNA on chromosomes.pdf
pdf("Location of eccDNA on chromosomes.pdf",width=16,height = 12)
kp <- plotKaryotype()
kpPlotDensity(kp, data = gr_example1,col="blue",r0 = -0.5,r1 = -1.2)
kpPlotDensity(kp, data = gr_example2,col="red",r0 = 0,r1 = 0.7)
dev.off()

```

 2. 结果演示

![eccDNA在染色体上的分布](https://cdn.nlark.com/yuque/0/2023/png/38459666/1698247491577-2fd0ad50-7fcb-48d9-b347-b71301cffdcb.png)

eccDNA在所有染色体上均有分布，但分布并不均匀

### <font style="color:rgb(79, 79, 79);">方法2：使用RIdeogram 包</font>
教程链接：[https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html](https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html)

这个包的优点是安装方便，不需要BiocManager，并且依赖其他包的数量也比较少



1. 获得绘图所需要的人类染色体长度及着丝粒(Centromere )

这个包自带的数据为hg38的人类染色体长度及着丝粒(Centromere )，如果需要hg19的数据需要自行下载，下载，获取方法教程如下：

[https://www.jianshu.com/p/208e479bc36a](https://www.jianshu.com/p/208e479bc36a)

获取后只保留染色体chr1~chr22，chrX和ChrY整理成如下的tab分隔符的txt文件

[human_karyotype.txt](https://www.yuque.com/attachments/yuque/0/2023/txt/38459666/1699501987441-032787f9-70bd-40c1-a6e4-5bc032de0b6a.txt)



2. 使用linux 的awk将bed转为gff3

```bash
awk -v OFS="\t" '{print $1,"HAVANA","gene",$2+1,$3,".","+",".","."}' 0h.bed > 0h.gff3
awk -v OFS="\t" '{print $1,"HAVANA","gene",$2+1,$3,".","+",".","."}' 16h.bed > 16h.gff3
```

3. 绘图

```bash
# 设置工作路径
setwd(dir = "D:\\HepG2\\代码\\srtp")

# 安装所需要的包
install.packages("RIdeogram")

# 导入所需数据
human_karyotype <- read.csv("human_karyotype.txt",sep = "\t",header = T)

# 以1000000bp为窗口，计算每个窗口内的eccDNA数量
Control_density <- GFFex(input = "0h.gff3", karyotype = "human_karyotype.txt", feature = "gene", window = 1000000)
DOX_density <- GFFex(input = "16h.gff3", karyotype = "human_karyotype.txt", feature = "gene", window = 1000000)

# 使用ideogram包绘图
ideogram(karyotype = human_karyotype, overlaid = Control_density, label = DOX_density, label_type = "heatmap", colorset1 = c("#f7f7f7", "#2c7fb8"), colorset2 = c("#f7f7f7", "#e34a33"))
convertSVG("chromosome.svg", device = "png")
```

4. 结果演示

![](https://cdn.nlark.com/yuque/0/2023/png/38459666/1699502865327-1514b99f-8027-4cf7-951e-fd41f8539463.png)

## bigwig（或称bw）文件的使用方法
bigwig数据格式介绍可访问下面的blog，相比于bam和wig文件，这种格式的存储所需空间小，软件处理较快。360云盘已经上传bigwig文件，这一章的目的主要是介绍bigwig是怎么来的，以及后续的使用方法。

[生信格式 | bigwig，bw （基因组浏览器绘制）_bigwig文件-CSDN博客](https://blog.csdn.net/u011262253/article/details/109367884)

### 如何获得bigwig文件
下面主要介绍使用deeptools中的bamCoverage将bam文件转为bigwig，wig或bedgraph转为bigwig可以参照上面的blog

### 使用deeptools的bamCoverage或bamCompare工具
```bash
# 安装deeptools
## deeptools的安装方法同样可以使用conda
conda install -c bioconda deeptools
## 如果失败可以新建一个环境，然后在新环境中安装
conda create -n deeptools
conda activate deeptools
conda install -c bioconda deeptools

# bam文件转为bigwig
## 需要构建bam的索引
samtools index test.bam

# 下面代码参数的简介
## -p设置线程数（使用‘cat /proc/cpuinfo |grep 'processor'|sort -u|wc -l’查看线程数）
## --normalizeUsing CPM，均一化数据格式的方式为“CPM”，还可以使用RPKM,BPM,RPGC,None
## --binSize 50 （软件的默认值就是50，所以不加也是可以的），计算数据使用的窗口大小。
## 可以使用bamCoverage -h查看这些参数的含义。
## "\"为脚本的换行符，可以使代码更清楚。不加"\"相当
## 于"bamCoverage -p 8 --normalizeUsing CPM --bin#Size 50 -b test.bam -o test.bw"

bamCoverage -p 8 \
--normalizeUsing CPM \
--binSize 50 \
-b test.bam \
-o test.bw

# 如果是有IP组和input组，如Chip-seq可以使用bamCompare
## Chip-seq介绍：https://zhuanlan.zhihu.com/p/279354841
## 同样先构建bam的索引，如果已经构建可以忽略
samtools index test_IP.bam
samtools index test_input.bam

bamCompare -p 8 \
-b1 test_IP.bam \
-b2 test_input.bam \
--binSize 50 \
--normalizeUsing CPM \
-o test.bw
```

### bigwig的使用方法
1. 使用deeptools的multiBigwigSummary工具对特定bed区域定量

```bash
# 参数介绍
## 设置BED-file对特定的bed文件进行定量，下面输入的example1.bed有4307行，则multiBigwigSummary会根据输入的bigwig文件对每一行的染色体区间进行定量
## --numberOfProcessors设置线程数
## --bwfiles输入bigwig文件
## --BED输入对应的bed文件
## --labels设置对应组学的显示标签
## scores_example1_ecc.tab为一个tab分隔符的文件，每一行对应提供bed文件的每一行
multiBigwigSummary BED-file \
--numberOfProcessors 8 \
--bwfiles test1.bw test2.bw test3.bw \
--BED example1.bed \
--labels H3K27ac ATAC RNA-seq \
-out scores_example1_ecc.npz --outRawCounts scores_example1_ecc.tab
```

下面是得到的示例数据，包含example1.bed和example2.bed的计算结果，可以使用excel打开，可以使用这些值进行统计分析。

[scores_ecc.zip](https://www.yuque.com/attachments/yuque/0/2023/zip/38459666/1699192547628-a2911840-ecbf-4608-952b-25b8331ec5a4.zip)



2. 使用computeMatrix和plotProfile可视化特定bed区域的数值

```bash
# 使用deeptools的computeMatrix的工具
## --beforeRegionStartLength,bed区间的起点向上游延伸的距离
## --regionBodyLength将bed区间缩放到1000bp
## --afterRegionStartLength，bed区间的终点向下游延伸的距离
## --numberOfProcessors 8设置线程数量，数量越高运行速度越快
## 计算数值的窗口，窗口越大绘图的结果约平滑，默认为10bp
computeMatrix scale-regions -S test1.bw \
-R example1.bed example2.bed \
--beforeRegionStartLength 0 \
--regionBodyLength 1000 \
--afterRegionStartLength 0 \
--numberOfProcessors 8 \
--binSize 10 \
-o test1.tab.gz

plotProfile \
-m test1.tab.gz \
-out test1.pdf \
--startLabel "ecc_start" \
--endLabel "ecc_end" \
--plotFileFormat "pdf"
```

结果示例如下

![](https://cdn.nlark.com/yuque/0/2023/png/38459666/1699195838603-b654f023-7d34-48a8-8868-2d6012d7c6eb.png)



## 


## Chip-seq call peak
Chip-seq的call peak过程和ATAC-seq的过程类似，除了在macs2 call peak的过程中，并没有使用参数"--nomodel --shift -100 --extsize 200"其余是一样的，可以在macs2_callpeak文件夹下找到以".narrowPeak"为结尾的文件，这个就是我们需要的peak区域，这些Chip-seq H3k27ac的peak区域意味着在基因组中有很高乙酰化水平的区域。同样的，文件已经上传到360云盘。

```bash
#进入工作路径
cd /newdisk/yuguo/Chip-seq_hepg2

#文件路径
data_dir=/newdisk/Sequencing_data/Dox_hepg2_p53/Chip_seq/old

#文件路径
data_dir=/newdisk/Sequencing_data/Dox_hepg2_p53/Chip_seq/old

#创建文件夹存放不同的文件
mkdir -p macs2_callpeak
mkdir -p bwfile
mkdir -p bamfile
mkdir -p txtfile

#创建文件软连接，方便后续对文件操作
ls $data_dir|while read id;do ln -s $data_dir/$id $id;done
 
#参考基因组路径
genome_index=~/genome_index/bowtie2/hg19

#创建文件名txt
ls *_1.fq.gz|while read id;do echo $( basename $id _1.fq.gz ) >> txtfile/file_name.txt;done
cat txtfile/file_name.txt | grep IP > 1
cat txtfile/file_name.txt | grep input > 2
paste 1 2 > txtfile/IP_and_input.txt
rm 1 2

#激活环境
source activate atac 

cat txtfile/file_name.txt | while read id
do
## cat txtfile/file_name.txt
## Input_0h_clean
## Input_16h_clean
## IP_0h_clean
## IP_16h_clean
#因为测序公司质控过了，所以跳过质控这一步

#比对到参考基因组
bowtie2 -x ${genome_index} -1 ${id}_1.fq.gz  -2 ${id}_2.fq.gz -p 15 | samtools view /dev/stdin -@ 15 -u |  samtools sort /dev/stdin -@ 15 -o bamfile/${id}.sorted.bam

#去除低质量序列
#-q 30 去除质量低于30的reads -f 2 保留正确比对到参考基因组的reads对
samtools view bamfile/${id}.sorted.bam -h -f 2 -q 30 -@ 10 -b > bamfile/${id}.sorted.filt.bam

#picard标记重复
picard MarkDuplicates I=bamfile/${id}.sorted.filt.bam O=bamfile/${id}.sorted.dedup.bam M=txtfile/${id}.dup_metrics.txt REMOVE_DUPLICATES=true

#索引bam文件为后续创建bw文件做准备
samtools index bamfile/${id}.sorted.dedup.bam

#用bedtools将bam文件转换为bedpe文件
samtools sort -n bamfile/${id}.sorted.dedup.bam -@ 12|bedtools bamtobed -bedpe -mate1 -i /dev/stdin | gzip -nc > bamfile/${id}.bedpe.gz

#创建tagAlign文件
zcat bamfile/${id}.bedpe.gz | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > bamfile/${id}.tagAlign.gz

done

cat txtfile/IP_and_input.txt | while read id
## cat txtfile/IP_and_input.txt
## IP_0h_clean	Input_0h_clean
## IP_16h_clean	Input_16h_clean
do
var=( $id )
treat=${var[0]}
control=${var[1]}
prefix=$( basename $treat _IP_clean )

#macs2 callpeak
macs2 callpeak -t bamfile/${treat}.tagAlign.gz -c bamfile/${control}.tagAlign.gz --format BED -g hs -q 0.05 -n ${prefix} --keep-dup all --SPMR --outdir macs2_callpeak

#创建可视化轨道
bamCompare -p 10 -b1 bamfile/${treat}.sorted.dedup.bam \
-b2 bamfile/${control}.sorted.dedup.bam \
--binSize 50 \
--normalizeUsing CPM \
-o $prefix.bw

done
```



## RNA-seq 数据处理
### 序列比对
```bash
conda install -c bioconda hisat2

# 构建索引
# hisat2-build genome.fa genome.fa

# 下载索引
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/hg19.tar.gz
tar -zxvf *.tar.gz
rm -rf *.tar.gz

##shell脚本；使用for循环嵌套
for o in 0 8 16
do

for i in 1 2 3
do
##-p 线程； -t 时间； --dta 转录组； -x 参考基因组路径； -1 双端测序结果的第一个文件； -2 双端测序结果的第二个文件； -S 输出.sam文件
hisat2 -p 10 \
-t --dta \ 
-x ~/genome_index/hisat_hg19/hg19 \
-1 /newdisk/yuguo/RNA-seq_hepg2/data/fastq/HepG2___用Dox处理"$o"h-"$i"_1.clean.fq.gz \
-2 /newdisk/yuguo/RNA-seq_hepg2/data/fastq/HepG2___用Dox处理"$o"h-"$i"_2.clean.fq.gz \
-S HepG2___用Dox处理"$o"h-"$i".sam

##第一步view将比对后的sam文件转换成bam文件。-@ 线程； -S 后面跟的是sam文件的路径；-b 指定输出的文件为bam，后面跟输出的路径；最后重定向写入bam文件
##第二步sort将所有的bam文件按默认的染色体位置进行排序。-o是输出文件名
##第三步index将所有的排序文件建立索引，生成的索引文件是以bai为后缀的
samtools view -@ 10 -S "$o"h-"$i".sam -b > "$o"h-"$i".bam
samtools sort "$o"h-"$i".bam -o "$o"h-"$i"_sorted.bam
samtools index /newdisk/yuguo/RNA-seq_hepg2/data/"$o"h-"$i"_sorted.bam

done
done
```



### 转录本拼接
转录本拼接用于对全基因组范围内的所有基因进行定量，结果以.tab结尾，已经上传到360云盘"`我的文件\2022年阿霉素处理hepg2细胞测序数据\Circle-Seq\RNA-seq数据`"文件夹内，"stringtieTPM.xlsx"是经过整理后的结果

```bash
# 安装stringtie软件
conda search -c bioconda stringtie

# 下载所需的基因注释文件
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.basic.annotation.gff3.gz
gunzip gencode.v42.basic.annotation.gff3.gz

## 
stringtie  -p 10 "$o"h-"$i"_sorted.bam -A "$o"h-"$i".tab -C "$o"h-"$i"_cov_refs.gtf -G /home/yuguo/annotation/gencode.v42lift37.basic.annotation.gff3 -o "$o"h-"$i".gtf

```



### reads定量
因为常见的差异分析软件会使用原始的reads数量作为输入，所以进行基因的差异分析时需要使用reads数量对基因进行定量，常用的软件为featureCounts。结果已经上传至360云盘"`我的文件\2022年阿霉素处理hepg2细胞测序数据\Circle-Seq\RNA-seq数据`"，命名为featurecount结果.txt。

<font style="color:rgb(51, 51, 51);">featurecounts结果解读：</font>

<font style="color:rgb(51, 51, 51);">Geneid 代表统计的meta-features 的名称（这里是该基因在ensemble数据库中的id），Chr(染色体编号)，Start，End 染色体上的位置，Strand 正负链，Length 所有exon区间长度的总和，最后6列是输入文件的名称，代表的是这个基因 的count 值，即表达量。</font>

<font style="color:rgb(51, 51, 51);">若第一个geneid在gtf 文件中有6个exon（外显子）的记录，就是说有6个features , 即对应6个Chr, Start, End, Strand，用分号分隔6个值。</font>

```bash
featureCounts -T 10 -p --countReadPairs -t exon -g gene_id -a ~/annotation/Homo_sapiens.GRCh37.87.gtf -o all.id.txt *_sorted.bam
##-T 线程数； -p --countReadPairs 如果指定，会统计fragment而不统计read； -t 设置feature-type，-t必须是gtf中有的feature，同时read只有落到这些feature上才会被统计到，默认是“exon”； -g 必须是gtf中提供的id identifier，默认为gene_id； -a 参考注释文件的路径； -o 输出文件名
```



### 基因的差异分析
<font style="color:rgb(51, 51, 51);">把featurecount结果.txt文件的多余数据删掉（只留下geneid和read counts），删除重复数据（Excel--数据--删除重复值），保存为0-16h.csv。差异分析的结果已经上传到360云盘"</font>`<font style="color:rgb(51, 51, 51);">我的文件\2022年阿霉素处理hepg2细胞测序数据\Circle-Seq\RNA-seq数据</font>`<font style="color:rgb(51, 51, 51);">"文件夹内，"</font>`<font style="color:rgb(51, 51, 51);">0h-16h-id.csv</font>`<font style="color:rgb(51, 51, 51);">"为注释过且padj<0.05的具有显著差异的基因信息,"</font>`<font style="color:rgb(51, 51, 51);">0h-16h_results.csv</font>`<font style="color:rgb(51, 51, 51);">"为所有基因的差异分析结果。</font>

```bash
##DESeq2差异分析和gene_id的注释信息
BiocManager::install(c("R.utils","DESeq2","tidyverse","biomaRt","curl"))
setwd("D:\\HepG2\\代码\\srtp")
library("R.utils")
library("DESeq2")
library("tidyverse")
mycounts<-read.csv("0-16.csv")

##第一列geneid设为行名
rownames(mycounts)<-mycounts[,1]

##删除第一列的内容
mycounts<-mycounts[,-1]

##factor()因子；rep()重复函数，重复3次
condition <- factor(c(rep("0h",3),rep("16h",3)), levels = c("0h","16h"))

##将mycounts中的列名作为数据框的行名，并重命名为0h、0h...
colData <- data.frame(row.names=colnames(mycounts), condition)

##检查是否有NA(缺失值)，0不是缺失值
mycounts[is.na(mycounts)] <- 0

##构建dds矩阵，~为构建公式，左边为因变量，右边为自变量
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds) ##标准化

##差异分析结果
res = results(dds, contrast=c("condition", "16h", "0h")) 

##按pvalue值从小到大排序
res = res[order(res$pvalue),]
write.csv(res,file="0h-16h_results.csv")

##table()函数提取padj<0.05的值，进行基因注释，即生成gene_id文件
table(res$padj<0.05)

##subset()从res中筛选数据；abs()绝对值
diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

##dim()函数输出一个矩阵的行数，紧跟着输出这个矩阵的列数
dim(diff_gene_deseq2)

library('biomaRt')
library("curl") 

##启用biomaRt包中的数据库
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(diff_gene_deseq2)

##getBM()函数各种ID转换，并获得其注释信息
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                    filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
ensembl_gene_id<-rownames(diff_gene_deseq2)

## cbind()函数将多个向量合成一个数据框
diff_gene_deseq2<-cbind(ensembl_gene_id,diff_gene_deseq2)
colnames(diff_gene_deseq2)[1]<-c("ensembl_gene_id")

##merge()函数合并数据框，以geneid合并
diff_name<-merge(diff_gene_deseq2,mms_symbols,by="ensembl_gene_id")
write.csv(diff_name,file="0h-16h-id.csv")
```

## 断点核苷酸可视化
eccDNA的断点可能具有核苷酸特征，如下图（start和end之间的区域为eccDNA区域，start和end即为断点，断点首尾连接成环后形成eccDNA）。

![](https://cdn.nlark.com/yuque/0/2023/png/38459666/1700302753761-01273635-6244-440d-87fc-c05105e96ab6.png)

> <font style="color:rgb(33, 33, 33);">Sin STK, Jiang P, Deng J, Ji L, Cheng SH, Dutta A, Leung TY, Chan KCA, Chiu RWK, Lo YMD. Identification and characterization of extrachromosomal circular DNA in maternal plasma. Proc Natl Acad Sci U S A. 2020 Jan 21;117(3):1658-1665. doi: 10.1073/pnas.1914949117. Epub 2020 Jan 3. PMID: 31900366; PMCID: PMC6983429.</font>
>

**<font style="color:rgb(33, 33, 33);">绘图思路：</font>**

1. <font style="color:rgb(33, 33, 33);">获得eccDNA断点上下游10bp的bed文件</font>
2. <font style="color:rgb(33, 33, 33);">使用bedtools getfasta获得断点上下游10bp的核苷酸序列</font>
3. <font style="color:rgb(33, 33, 33);">使用R包ggseqlogo绘图</font>



### 需要的准备文件或软件
1. <font style="color:rgb(33, 33, 33);">示例bed文件</font>

[0h.zip](https://www.yuque.com/attachments/yuque/0/2023/zip/38459666/1700303267309-114b60de-732a-4d69-85ac-04add11bca71.zip)

2. 参考基因组的长度，下载方式如下(这里是ucsc版本）：

[如何获取染色体长度](https://zhuanlan.zhihu.com/p/343564985#:~:text=UCSC%E4%B8%BB%E9%A1%B5%E4%B8%8A%E7%82%B9%E5%87%BBDownloads-%3EGenome%20Data%202.%E9%80%89%E6%8B%A9%E7%9B%B8%E5%BA%94%E7%9A%84%E7%89%A9%E7%A7%8D%EF%BC%8C%E8%BF%99%E9%87%8C%E9%80%89%E6%8B%A9human%203.%E7%82%B9%E5%87%BBGenome%20sequence%20files%20and%20select,genome%20sequence%20files%20and%20select%20annotations%204.%20%E4%B8%8B%E8%BD%BDhg38.chrom.sizes%E8%BF%99%E4%B8%AA%E6%96%87%E4%BB%B6%E5%8D%B3%E5%8F%AF%EF%BC%8C%E8%B7%9F%E4%B8%8A%E9%9D%A2%E4%B8%80%E7%A7%8D%E6%96%B9%E6%B3%95%E5%BE%97%E5%88%B0%E7%9A%84%E7%BB%93%E6%9E%9C%E6%98%AF%E4%B8%80%E6%A0%B7%E7%9A%84%E3%80%82)

或命令行：`wget [https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes)`

[hg19.chrom.zip](https://www.yuque.com/attachments/yuque/0/2023/zip/38459666/1700304171962-d08e6061-b78a-49cf-ac19-8831b86b8733.zip)

3. 参考基因组fasta序列(这里是ucsc版本）

`wget -c [https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)`



4. bedtools软件

```bash
# 使用conda安装
conda install -c bioconda bedtools

# 或直接下载编译后的bedtools软件
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make
```



5. samtools 软件

```bash
# 使用conda安装
conda install -c bioconda samtools
```

### 代码
1. 先获得<font style="color:rgb(33, 33, 33);">eccDNA断点上下游10bp的bed文件</font>

```bash
# 获得参考基因组的长度
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

# 获得start扩展10bp后的文件
##  "awk -v OFS="\t" '$3-$2>20 {print $1,$2-10,$2+10}' 0h.bed"，删除长度小于20bp的eccDNA，并将起始位点（start）向左和向右移动10bp。如：一个eccDNA的起始点为chr1	50	100，则生成的上游文件为：chr1	40	60。20bp的删除是为了更严谨一些，当然也可以不删除。
## "bedtools slop -i /dev/stdin -g hg19.chrom.sizes -b 0" 将移动后的位置修正在染色体的范围内,“-b 0” 将bed的start向左移动0bp，将end向右移动0bp。如：有一个eccDNA的长度是chr1	5	30，起点向左向右移动10bp后为chr1 -5	15，起点小于0，这明显不合理，修正后为：chr1	0	15
 awk -v OFS="\t" '$3-$2>20 {print $1,$2-10,$2+10}' 0h.bed | bedtools slop -i /dev/stdin -g hg19.chrom.sizes -b 0 > start_flanking_10bp.bed

 # 获得end扩展10bp后的文件
 awk -v OFS="\t" '$3-$2>20 {print $1,$3-10,$3+10}' 0h.bed | bedtools slop -i /dev/stdin -g hg19.chrom.sizes -b 0 > end_flanking_10bp.bed
```



2. <font style="color:rgb(33, 33, 33);">使用bedtools getfasta获得断点上下游10bp的核苷酸序列</font>

```bash
# 解压hg19.fa.gz,-d解压后将删除hg19.fa.gz
gzip -d hg19.fa.gz

# 对hg19.fa构建索引(会产生一个"hg19.fa.fai"文件）
samtools faidx hg19.fa

# 获得对应的核苷酸序列，-fi参考基因组fasta序列，-bed输入的bed文件，-tab生成将以tab分隔而不是FASTA格式，-fo输出文件
bedtools getfasta -fi hg19.fa -bed start_flanking_10bp.bed -tab -fo start_flanking_10bp.txt
bedtools getfasta -fi hg19.fa -bed end_flanking_10bp.bed -tab -fo end_flanking_10bp.txt

# head start_flanking_10bp.txt
chr1:982254-982274      GGCTTCCTGCTACAACTCCG
chr1:1867378-1867398    aggtgtagtggctcatgcct
```

[flanking_10bp.zip](https://www.yuque.com/attachments/yuque/0/2023/zip/38459666/1700308593946-cf834e5e-8630-4092-a8d7-ada04211b420.zip)

3. 绘制图形

[R包ggseqlogo |绘制序列分析图](https://www.jianshu.com/p/541293447b7a)

```bash
# 设置工作路径
setwd(dir = "D:\\HepG2\\代码\\srtp")

# 安装ggseqlogo
#直接从CRAN中安装
install.packages("ggseqlogo")
#从GitHub中安装
devtools::install.github("omarwagih/ggseqlogo")

#加载包
library(ggplot2)
library(ggseqlogo)

#加载数据
start <- read.table("start_flanking_10bp.txt", header=F, row.names=NULL)
end <- read.table("end_flanking_10bp.txt", header=F, row.names=NULL)

> head(start)
                    V1                   V2
1   chr1:982254-982274 GGCTTCCTGCTACAACTCCG
2 chr1:1867378-1867398 aggtgtagtggctcatgcct
3 chr1:2203702-2203722 TAAATGTCTTAGTAGCTTAA
4 chr1:2226450-2226470 GTCAGGCCCCTGGCTCAGAT
5 chr1:3199335-3199355 CCGTCTGAAGCCTGGCGAAC
6 chr1:3279159-3279179 TGCCGAGGGCCTTCGGGTCA

# 绘图
ggplot()+geom_logo(start$V2)+theme_logo()
ggplot()+geom_logo(end$V2)+theme_logo()

# 更改主题/保存
library(ggprism)
theme=theme_prism(base_size = 14,base_line_size = 1,base_fontface = "plain")

ggplot()+geom_logo(start$V2)+theme_logo()+theme
ggsave(filename = "start_flanking_10bp.pdf",device = "pdf",width = 8,height = 4)

ggplot()+geom_logo(e$V2)+theme_logo()+theme
ggsave(filename = "end_flanking_10bp.pdf",device = "pdf",width = 8,height = 4)
```

### 结果展示
![画板](https://cdn.nlark.com/yuque/0/2023/jpeg/38459666/1700310457913-07371dd0-a103-48b1-bea6-91dd724354fc.jpeg)



### 序列标识图释义
[3.12寻找保守区域-01-序列标识图WebLogo.pdf](https://www.yuque.com/attachments/yuque/0/2023/pdf/38459666/1701830641499-0069102e-1101-43df-a861-068c3eaff4cf.pdf)

## eccDNA的数量与其他因素的关联
文献显示eccDNA出现在GC值高，基因密度高的区域（<font style="color:rgb(33, 33, 33);">PMID: 26051933</font>）。下方是原文的图，图中显示了eccDNA的数量在染色体10的分布，并标注了GC值和基因密度，可以从图中看到eccDNA的数量和GC值确实在多个组织中显现出较高的一致性。可及将此图扩展至多个组学特征和eccDNA数量的关联性。

![](https://cdn.nlark.com/yuque/0/2023/png/38459666/1700448116257-d7476a1c-6847-4d17-9dd2-4cab53a4ed47.png)

**绘图思路**

1. 获得染色体每Mb的eccDNA数量以及其他数据的值（基因数量，GC含量，甲基化含量，H3K27ac含量，ATAC-seq含量，RNA-seq含量）
2. 将每mb的eccDNA数量与其余数据进行相关性分析，并可视化
3. 挑选染色体进行可视化展示相关性。



### 需要的准备文件或软件
代码运行在R或linux平台

1. eccDNA

[ecc_hg38.bed.zip](https://www.yuque.com/attachments/yuque/0/2023/zip/38459666/1700448999903-352184b2-1012-467c-8078-665c327ae7d1.zip)（包含两个eccDNA文件，hg38版本）

2. RIdeogram包

可以使用此包的GFFex功能获得每Mb的eccDNA的数量和基因数量

> 教程：[https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html](https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html)
>

3. GC数据

下载GC含量的bigwig文件（bin值为5）

`wget -c [https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.gc5Base.bw](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.gc5Base.bw)`

4. 甲基化数据

# 下载负链bw文件

`wget -c [https://www.encodeproject.org/files/ENCFF464CWC/@@download/ENCFF464CWC.bigWig](https://www.encodeproject.org/files/ENCFF464CWC/@@download/ENCFF464CWC.bigWig)`

# 下载正链bw文件

`wget -c [https://www.encodeproject.org/files/ENCFF884UAE/@@download/ENCFF884UAE.bigWig](https://www.encodeproject.org/files/ENCFF884UAE/@@download/ENCFF884UAE.bigWig)`

5. H3K27ac，ATAC-seq，RNA-seq（hg38版本）

已上传至360云盘，名称为

> "0h_ATAC_hg38.bw"，"0h_H3K27ac_hg38.bw"，"0h_mRNA_hg38.bw"，"16h_ATAC_hg38.bw"，"16h_H3K27ac_hg38.bw"，"16h_mRNA_hg38.bw"
>

6. deeptools

`conda install -c bioconda deeptools`

7. 基因注释文件

`wget -c [https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gff3.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gff3.gz)`

[gencode.v44.basic.annotation.gff3.gz](https://www.yuque.com/attachments/yuque/0/2023/gz/38459666/1700458575659-c7a90192-55fa-4a29-9e04-a5a76a6f904b.gz)

8. python

[全网最详细的Python安装教程（Windows）](https://zhuanlan.zhihu.com/p/344887837)

### 代码
1. **获得染色体每1Mb的eccDNA数量以及其他数据的值**

**使用linux的awk工具将bed转为gff3文件**

```bash
awk -v OFS="\t" '{print $1,"HAVANA","gene",$2+1,$3,".","+",".","."}' 0h.hg38.bed > 0h.hg38.gff3
awk -v OFS="\t" '{print $1,"HAVANA","gene",$2+1,$3,".","+",".","."}' 16h.hg38.bed > 16h.hg38.gff3
```



**获得每Mb eccDNA的数量和基因的数量，使用R包RIdeogram**

```r
# setwd()

library(RIdeogram)
# 该包自带了人类染色体和着丝粒的位置
data(human_karyotype, package="RIdeogram")

# 注意该包的染色体是"1","2","3"类型，之后要转化为"chr1","chr2","chr3"类型
> head(human_karyotype)
  Chr Start       End  CE_start    CE_end
1   1     0 248956422 122026459 124932724
2   2     0 242193529  92188145  94090557
3   3     0 198295559  90772458  93655574
4   4     0 190214555  49712061  51743951
5   5     0 181538259  46485900  50059807
6   6     0 170805979  58553888  59829934

# 转化为"chr1","chr2","chr3"类型
human_karyotype$Chr <- paste("chr",human_karyotype$Chr,sep = "")
> head(human_karyotype)
   Chr Start       End  CE_start    CE_end
1 chr1     0 248956422 122026459 124932724
2 chr2     0 242193529  92188145  94090557
3 chr3     0 198295559  90772458  93655574
4 chr4     0 190214555  49712061  51743951
5 chr5     0 181538259  46485900  50059807
6 chr6     0 170805979  58553888  59829934

将其保存一份，之后能用到
write.table(x = human_karyotype,file = "human_karyotype_hg38.txt",quote = F,sep = "\t",row.names = FALSE)

# 以1mb为单位计算数量
Control_density <- GFFex(input =  "0h.hg38.gff3", karyotype = "human_karyotype_hg38.txt", feature = "gene", window = 1000000)
DOX_density <- GFFex(input = "16h.hg38.gff3", karyotype = "human_karyotype_hg38.txt", feature = "gene", window = 1000000)
gene_gensity <- GFFex(input = "gencode.v44.basic.annotation.gff3.gz", karyotype = "human_karyotype_hg38.txt", feature = "gene", window = 1000000)

#取消科学计数法
options(scipen = 200)

# 保存文件
write.table(x = Control_density,file = "Control_density_1mb.txt",quote = F,sep = "\t",row.names = FALSE,col.names = F)
write.table(x = DOX_density,file = "DOX_density_1mb.txt",quote = F,sep = "\t",row.names = FALSE,col.names = F)
write.table(x = gene_gensity,file = "gene_gensity_1mb.txt",quote = F,sep = "\t",row.names = FALSE,col.names = F)

```



2. **将1mb的eccDNA数量与其余数据进行相关性分析，并通过linux软件deeptools可视化**

```bash
# 文件的路径
mrna_0h="/mnt/d/HepG2/全局/hg38bw/0h_mRNA_hg38.bw"
mrna_16h="/mnt/d/HepG2/全局/hg38bw/16h_mRNA_hg38.bw"
atac_0h="/mnt/d/HepG2/全局/hg38bw/0h_ATAC_hg38.bw"
atac_16h="/mnt/d/HepG2/全局/hg38bw/16h_ATAC_hg38.bw"
h3k27ac_0h="/mnt/d/HepG2/全局/hg38bw/0h_H3K27ac_hg38.bw"
h3k27ac_16h="/mnt/d/HepG2/全局/hg38bw/16h_H3K27ac_hg38.bw"
GC="/mnt/d/HepG2/全局/hg38.gc5Base.bw"
meminus="/mnt/d/HepG2/全局/ENCFF464CWC.bigWig"
meplus="/mnt/d/HepG2/全局/ENCFF884UAE.bigWig"

# 计算
multiBigwigSummary BED-file \
--numberOfProcessors 10 \
 --bwfiles ${mrna_0h} ${atac_0h} ${h3k27ac_0h} ${mrna_16h} ${atac_16h} ${h3k27ac_16h} ${meplus} ${meminus} ${GC} \
 --BED Control_density_1mb.txt \
 --labels RNA-seq_0h ATAC-seq_0h H3K27ac_0h RNA-seq_16h ATAC-seq_16h H3K27ac_16h WGBS_plus_strand WGBS_minus_strand  WGBS_GC_content \
 -out scores_RAHMG_1mb.npz --outRawCounts scores_RAHMG_1mb.tab

```



**将"Control_density_1mb.txt""DOX_density_1mb.txt" "gene_gensity_1mb.txt"与scores_RAHMG_1mb.tab文件进行合并（excel）**

[all_data_1mb.txt](https://www.yuque.com/attachments/yuque/0/2023/txt/38459666/1700467403883-e1bdc4ba-ed70-4319-b97d-bd8c7b29911b.txt)



**将tab分隔的txt文件转为npz文件（python），用于之后的绘图**

```bash
# 安装numpy
# pip install numpy

import numpy as np

# 设置工作路径（也可以使用绝对路径）
# import os
# os.chdir("D:\HepG2\代码\srtp\")

# 读取tab分隔符文件
all= np.loadtxt(r"D:\HepG2\代码\srtp\all_data_1mb.txt")

# 将tab分隔的txt文件转为npz文件
np.savez(r"D:\HepG2\代码\srtp\all_data_1mb.npz",matrix=all,labels=np.array(['ecc_density_Control','ecc_density_DOX','gene_density','RNA-seq_Control','ATAC-seq_Control','H3K27ac_Control','RNA-seq_DOX','ATAC-seq_DOX','H3K27ac_DOX','WGBS_plus_strand','WGBS_minus_strand','WGBS_GC_content']))
```

[all_data_1mb.zip](https://www.yuque.com/attachments/yuque/0/2023/zip/38459666/1700467489324-3dd25c63-adfd-428b-804c-1e9816b1ab02.zip)



**使用Deeptools的plotCorrelation绘图**

```bash
# 热图形式
plotCorrelation \
    -in all_data_1mb.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation" \
    --whatToPlot heatmap --colorMap Reds --plotNumbers \
    --plotFileFormat "pdf" \
    -o heatmap_SpearmanCorr_all_data_1mb.pdf

# 或者是散点图
plotCorrelation \
    -in all_data_1mb.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation" \
    --whatToPlot scatterplot --colorMap Reds --plotNumbers \
    --plotFileFormat "pdf" \
    -o scatterplot_SpearmanCorr_all_data_1mb.pdf

    plotCorrelation \
    -in all_data_5mb.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation" \
    --whatToPlot heatmap --colorMap Reds --plotNumbers \
    --plotFileFormat "pdf" \
    -o heatmap_SpearmanCorr_all_data_5mb.pdf
```

[heatmap_SpearmanCorr_all_data_1mb.pdf](https://www.yuque.com/attachments/yuque/0/2023/pdf/38459666/1700468072807-fc5d9e80-db13-4c33-9859-a9baceffe862.pdf)

当富集到的eccDNA数量有限的时候使用1mb作为窗口，可能并不能很好的反映相关性，因为大部分窗口都是1个或0个eccDNA，这个时候可以提高窗口的值，比如说5mb。这个时候的窗口值可能更能反映实际的相关性情况

[heatmap_SpearmanCorr_all_data_5mb.pdf](https://www.yuque.com/attachments/yuque/0/2023/pdf/38459666/1700468630450-da94d244-b286-46b5-8111-42012e5f8630.pdf)



3. **挑选染色体进行可视化展示相关性**

这里选染色体20作为示例

[数据匹配1mb.csv](https://www.yuque.com/attachments/yuque/0/2023/csv/38459666/1700482149955-790b8102-a45b-4dd1-8123-9f73211967a1.csv)

```r
library(ggplot2)
library(cowplot)

data <- read.csv("数据匹配1mb.csv",header = T)
data <- subset(data,data$Chr=="chr20")

# 定义主题
theme=theme_cowplot()+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

# 折线图线的粗细
linewidth=1

# 绘制图像
p1 <- ggplot()+geom_line(data = data,aes(x=Start,y = ecc_density_DOX.),linewidth=linewidth)+labs(x = "")+theme
p2 <- ggplot()+geom_line(data = data,aes(x=Start,y = gene_density.),linewidth=linewidth)+labs(x = "")+theme
p3 <- ggplot()+geom_line(data = data,aes(x=Start,y = X.RNA.seq_16h.),linewidth=linewidth)+labs(x = "")+theme
p4 <- ggplot()+geom_line(data = data,aes(x=Start,y = X.ATAC.seq_16h.),linewidth=linewidth)+labs(x = "")+theme
p5 <- ggplot()+geom_line(data = data,aes(x=Start,y = X.H3K27ac_16h.),linewidth=linewidth)+labs(x = "")+theme
p6 <- ggplot()+geom_line(data = data,aes(x=Start,y = X.WGBS_plus_strand.),linewidth=linewidth)+labs(x = "")+theme
p7 <- ggplot()+geom_line(data = data,aes(x=Start,y = X.WGBS_minus_strand.),linewidth=linewidth)+labs(x = "")+theme
p8 <- ggplot()+geom_line(data = data,aes(x=Start,y = X.WGBS_GC_content.),linewidth=linewidth)+labs(x = "Chr20")+theme_cowplot()

# 合并图像，ncol=1 将所有图像放置为1列，align = "vh",以水平和垂直对齐图像
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,ncol=1,align = "vh")

# 保存为pdf
ggsave("Chr20_correlation.pdf",device = "pdf",width = 10,height = 8)
```



**绘图后经过简单修图后效果如下：**

[Chr20_correlation.pdf](https://www.yuque.com/attachments/yuque/0/2023/pdf/38459666/1700482184059-694df637-6834-44f2-94e5-41f695544acc.pdf)

## eccDNA注释
### 使用homer对eccDNA进行注释
eccDNA注释可以了解eccDNA区域内包含哪些基因和元件

工具：homer

安装软件需要的硬件信息：

+ <font style="color:rgb(0, 0, 0);">Unix-style operating system (UNIX/LINUX/Mac/Cygwin)</font>
+ <font style="color:rgb(0, 0, 0);">1 Gb of RAM (4+ Gb)</font>
+ <font style="color:rgb(0, 0, 0);">1 Gb of Hard Drive Space (>10Gb)</font>

****

**安装并配置homer**

```bash
# 使用conda安装homer
conda install -c bioconda homer

# 还需要配置homer软件，下载必要的注释文件。需要用到configureHomer.pl脚本，这这脚本在安装homer时也被下载了下来，下面的命令可以找到configureHomer.pl的位置
find -name configureHomer.pl

#运行结果如下：
##./miniconda3/share/homer/configureHomer.pl

# 下面使用脚本下载参考基因组的注释文件，我下载的是hg19，可以按需要下载hg38，执行脚本需要一些时间
perl ./miniconda3/share/homer/configureHomer.pl -install hg19
```

****

**对eccDNA进行注释**

```bash
# 使用annotatePeaks脚本注释，--Stat example1.anno.stat对结果进行统计。结果有两个文件:example1.anno.stat和example1.anno.xls
annotatePeaks.pl example1.bed hg19 --Stat example1.anno.stat > example1.anno.xls
```

****

**可视化注释结果**

可以使用R或prism软件，推荐prism软件，可以自行搜索使用方法，绘图示例：

![](https://cdn.nlark.com/yuque/0/2023/png/38459666/1699187599078-585dc72e-633f-48a7-923c-3d969b7805ef.png)

蓝色example1，红色example2



### deeptools可视化eccDNA在基因区域的分布
**下载基因注释文件hg19版本**

[Data import into Galaxy — deepTools 3.5.4 documentation](https://deeptools.readthedocs.io/en/develop/content/help_galaxy_dataup.html?highlight=strand#download-annotation-files-from-public-data-bases)



**提取基因信息**

```bash
zcat gencode.v44lift37.basic.annotation.gff3.gz | awk -v OFS="\t" '$3=="gene" {print}' | sed 's/;/\t/g' | awk -F '\t' -v OFS="\t" '{gsub(/ID=/,"");gsub(/gene_type=/,"");print $1,$4,$5,$9,$11,$7}' > hg19_gene.bed
```

****

**将bed文件转为bw文件**

[ecc_hg19.bed.zip](https://www.yuque.com/attachments/yuque/0/2023/zip/38459666/1700572100555-2c971620-47c2-4154-8849-02aa01cf06f9.zip)

```bash
for name in 0h 16h
do
# bed3转为bed6
awk -v OFS="\t" '{print $1,$2,$3,".",".","."}' ${name}.bed > ${name}.bed6
# bed转为bam
# -i 输入文件 -g 基因组每条染色体的长度 -O BAM 输出格式 -o 输出文件的名称
bedToBam -i ${name}.bed6 -g hg19.chrom.sizes | samtools sort -O BAM -o ${name}.bam
# 构建bam索引
samtools index ${name}.bam
# bam转为bw文件
# -b 输入 -o 输出
bamCoverage -b ${name}.bam -o ${name}.bw --numberOfProcessors 10
done

```

****

**使用deeptools可视化**

```bash
# -S 输入的big文件 -R 输入的BED/GTF文件 --beforeRegionStartLength 提供的BED/GTF的起始点向上游延展的长度 --regionBodyLength 将提供的BED/GTF的长度统一缩放至/延展至的长度 --afterRegionStartLength 提供的BED/GTF的终点向下游延展的长度--skipZeros 舍去计算分数为零的BED/GTF区域
computeMatrix scale-regions -S 0h.bw 16h.bw -R ~/annotation/RefSeq.bed --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros --numberOfProcessors 10 -o ecc.tab.gz

# 绘制折线图
# -m 输入 --perGroup 以BED/GTF为分组将所有bw绘制到每个组中，若不添加则会以为bw分组将所有BED/GTF绘制到每个组中
plotProfile -m ecc.tab.gz --perGroup --plotTitle "Distribution of eccDNA on genes" --regionsLabel "" -o "Distribution of eccDNA on genes.pdf"
```



**结果展示**

结论：eccDNA在基因区域以及上下游3kb并没有出现偏好性分布

[Distribution of eccDNA on genes.pdf](https://www.yuque.com/attachments/yuque/0/2023/pdf/38459666/1700622995389-d76d4a4d-3781-4feb-aa6e-3fd0adb1821a.pdf)

### ChIPseeker可视化eccDNA在基因区域的分布
使用deeptools表现不出预期的结果，这里再换一种方法，看是否会出现预期的结果，这里我们使用R包ChIPseeker

```bash
# setwd("D:\\HepG2\\代码\\srtp")

# 安装ChIPseeker包
BiocManager::install("ChIPseeker")

# 加载包
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# 定义注释文件
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# 载入eccDNA的bed文件
peak <- readPeakFile("0h.bed")
peak2 <- readPeakFile("16h.bed")

# 将两个peak合为一个list
peak_list <- list(peak,peak2)

# 绘图并保存
pdf(file = "Distribution of eccDNA on genes_chipseeker.pdf",width = 8,height = 4)
plotPeakProf2(peak = peak_list, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 200,
              TxDb = txdb,ignore_strand = F,facet="row")
dev.off()
```

**结果展示**

可以看到，eccDNA的分布在基因的5'端减少了，这是我们想要的结果

[Distribution of eccDNA on gene.pdf](https://www.yuque.com/attachments/yuque/0/2023/pdf/38459666/1700655083757-024f9fbb-658f-40a4-8571-8a82faf50c59.pdf)



## eccDNA区域的模拟
若要生成基因组随机区域以及随机区域的fastq文件，可以使用Circle-Map

### 软件的安装和使用说明
```bash
# 使用conda新建一个环境
conda create -n circle-map

# 安装circle-map
conda install -c bioconda circle-map

# 激活虚拟环境
conda activate circle-map

# 查看Simulate参数的使用说明
Circle-Map Simulate -h
# -g 基因组fasta文件
# -N , --read-number 需要模拟的fastq序列的条目数量
# -o , --output 输出的文件名称
# -dir , --directory 输出的文件夹名称，默认为当前
# -b , --base-name fastq的名称前缀
# -s , --skip-region 跳过的区域
# -r , --read-length 模拟出的read的长度
# -m , --mean-insert  插入长度分布的平均值
# -c , --mean-coverage eccDNA区域的平均覆盖度
# -p , --processes 使用的线程数
# -v , --variants 如果为"true",引进突变
# -S , --substitution 对碱基进行替换的比例
# -I , --Indels 在基因组中引入 indels 的比例。默认值: 0.001
# -J , --java_memory Java 内存分配。默认值: -Xmx16g
# -e ,  --error 引入测序错误
# -i , --instrument 模拟Illumina测序仪器(Default HiSeq 2500)
# -ir , --insRate 首次读取插入率（默认值：0.00009）。需要 -e
# -ir2 , --insRate2 第二读取插入率（默认值：0.00015）。需要 -e
# -dr , --delRate 第一读取删除率（默认值：0.00011）。需要 -e
# -dr2 , --delRate2 第二读取删除率（默认值：0.00023）。需要 -e
```



### eccDNA的模拟过程
```bash
#！！！若要更改模拟eccDNA的大小可以仿造下面的操作方式
vim ~/anaconda3/envs/circle-map/lib/python3.6/site-packages/circlemap/simulations.py
# 使用vim编辑器中的查找命令搜索"circle_length = rd.randint"更改里面的"circle_length = rd.radint(150,350)",括号内的范围即为模拟的片段长度,默认值为150-350

genome_index=~/GRCh37_hg19/hg19.fa
genome_index=~/genome/hg38.fa
lengh=150
read_number=5000000
processes=8
# 想要屏蔽的区域，可以使测序的gap区域，可以到UCSC网站下载，还可以加上eccDNA区域
skip_region=~/annotation/gap.bed
bed_name="simulate.bed"
fa_name="simulate"
Circle-Map Simulate -g $genome_index -N $read_number -r $lengh -p $processes -s $skip_region -o $bed_name -b $fa_name

#hg19.chrom.sizes包含了常规染色体，我们提取包含这些染色体的前15000行
awk -v OFS="\t" 'NR==FNR{a[$1]; next} $1 in a' hg19.chrom.sizes simulate.bed | head -n 15000 > simulate.bed
```

## ATAC-seq 数据分析
### ATAC-seq call peak
> **Call Peak其实是一种计算方法，用于识别基因组中由测序得到的比对reads富集的区域。**  
其实通俗的来说Call peak就是把我们有蛋白富集到核酸时的区域给找出来。
>
> + 对于ChIPseq/Cut&Tag来说我们这个区域就是我们的转录因子所可以和DNA互作的区域；
> + 对于ATAC-seq来说这个就是特定条件下的开放染色质的区域；
> + 对于m6A测序来说就是发生m6A甲基化修饰的RNA区域。
>
> 链接：https://www.jianshu.com/p/d62b42ab0d51
>

可以在macs2_callpeak文件夹下找到以".narrowPeak"为结尾的文件，这个就是我们需要的peak区域，这些ATAC-seq的peak区域意味着在基因组中有很高染色质可及性的区域。文件已经上传到360云盘

```bash
## 创建conda atac所需要的conda环境
conda create -n atac

# 下载bowtie2软件用于序列比对，samtools用于对sam文件排序并转换为bam，picard用于区域比对序列中的重复，macs2用于计算peak区域的位置

conda install -c bioconda bowtie2
conda install -c bioconda samtools用于对sam文件
conda install -c bioconda picard
conda install -c bioconda macs2

## 构建bowtie2索引
bowtie2-build hg19.fa

## 使用等号定义变量，这里将"~/genome_index/bowtie2/hg19"赋值给genome_index，下面所有的genome_index都相当于"~/genome_index/bowtie2/hg19"
genome_index=~/genome_index/bowtie2/hg19
out=~/atac
cd $out

source activate atac

cat file_name.txt | while read id
do
## cat file_name.txt 
## HepG2_Dox_16h_Atac_clean
## HepG2_Dox_0h_Atac_clean
## "|"为管道符，用于将上一步得到数据传入下一步骤中，这里将输出的两行传给while函数，while函数对这两行逐一遍历，并将遍历的结果赋值给一个叫id的变量

#Alignment
## -x 输入的参考基因组
## -p 线程数量
## -u 在samtools中用于加速管道之间的数据传递
## -@ 5设置samtools的线程为5
bowtie2 -x ${genome_index} -1 ${id}_1.fq.gz  -2 ${id}_2.fq.gz -p 5 | samtools view -@ 5 -u -  |  samtools sort -@ 5 -o ${id}.sorted.bam

samtools index ${id}.sorted.bam
## -v 反选，-P"^chrM$"正则表达的模式
## xargs:将回车分隔变成空格分隔
## -h保留bam的头部信息
## -f 2 保留“read mapped in proper pair”
## 去除低质量序列（mapq>30）
## -b输出格式为bam
samtools idxstats ${id}.sorted.bam | cut -f 1 | grep -v -P "^chrM$" | xargs samtools view ${id}.sorted.bam -h -f 2 -q 30 -@ 5 -b> ${id}.sorted.rmChrM.bam

#去除重复
## O=""输出的bam文件名称
## M=""输出的矩阵的统计信息
picard MarkDuplicates I=${id}.sorted.rmChrM.bam O=${id}.sorted.rmChrM.dedup.bam M=${id}.dup_metrics.txt REMOVE_DUPLICATES=true

#用bedtools将bam文件转换为bedpe文件
samtools sort -n ${id}.sorted.rmChrM.dedup.bam -@ 5|bedtools bamtobed -bedpe -mate1 -i /dev/stdin | gzip -nc >${id}.bedpe.gz

#创建tagAlign文件
zcat ${id}.bedpe.gz | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > ${id}.tagAlign.gz

#macs2 callpeak
mkdir -p 
macs2 callpeak -t ${id}.tagAlign.gz --format BED -g hs --nomodel --shift -100 --extsize 200 -q 0.05 -n ${id} --keep-dup all --SPMR --outdir macs2_callpeak

samtools index ${id}.sorted.rmChrM.dedup.bam

bamCoverage -p 10 \
--normalizeUsing CPM \
--binSize 50 \
-b ${id}.sorted.rmChrM.dedup.bam \
-o ${id}.bw

done
```



### 基因区域的染色质可及性
使用ATAC-seq bigwig数据可以对基因区域的染色质可及性进行可视化

```bash
computeMatrix scale-regions -S 0h_atac.bw 16h_atac.bw -R ~/annotation/RefSeq.bed --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros --numberOfProcessors 10 -o ecc.tab.gz

# 绘制折现图
# -m 输入 --perGroup 以BED/GTF为分组将所有bw绘制到每个组中，若不添加则会以为bw分组将所有BED/GTF绘制到每个组中
plotProfile -m ecc.tab.gz --perGroup --plotTitle "Chromatin accessibility of gene regions"  --regionsLabel "" -o "Chromatin accessibility of gene regions.pdf" 
```

结果如下图：

![](https://cdn.nlark.com/yuque/0/2023/png/38459666/1701239944181-b0c5dbf7-8d01-44ff-bb30-7eeeefc74268.png)

## eccDNA区域的GC值可视化
使用deeptools的computeMatrix工具和plotProfile工具

```bash
# 下载GC含量的bigwig文件（bin值为5）

# hg38
# wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.gc5Base.bw

# hg19
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.gc5Base.bw

file="example.bed"
name=$(basename $file .bed)
computeMatrix scale-regions \
 -S hg19.gc5Base.bw \
 -R  ${name}.bed \
 -bs 10 \
--beforeRegionStartLength 1000 \
--regionBodyLength 1000 \
--afterRegionStartLength 1000 \
 --numberOfProcessors 8 \
 -out ${name}_gcflank1000.tab.gz

plotProfile \
 -m ${name}_gcflank1000.tab.gz \
 -out ${name}_gcflank1000.pdf \
--startLabel ecc_start \
--endLabel ecc_end \
--regionsLabel "${name}" \
--plotFileFormat "pdf"
```

