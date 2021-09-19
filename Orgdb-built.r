1.1 #读取数据
setwd("J:/C盘/Project/psc-langlab/task/apple/wangxiaofei/2021年8月19日-DREB2A/ko/ko")
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(AnnotationForge)

# 当你发现用read.table, read.csv出错时，用import是个不错的选择
egg <- rio::import('MM_a3a06bi4.emapper.annotations.tsv')
egg[egg==""] <- NA 

colnames(egg)

1.2 #提取基因ID
#重点是需要两列，一列是ID，一列是name（第二列不一定是真的name，但必须要有，比如看下面的示例）
gene_info <- egg %>%
  dplyr::select(GID = query, GENENAME = seed_ortholog) %>% na.omit()
  
1.3 #提取基因与GO的对应
gterms <- egg %>%
  dplyr::select(query, GOs) %>% na.omit()
  
gterms$GOs[which(gterms$GOs=="-")]=NA 
gterms <- na.omit(gterms)  
  
1.4 #再将gterms的每一行中GO按照逗号分开
# 之前的代码是用for循环，但是速度很感人，现在直接上sapply
library(stringr)
all_go_list=str_split(gterms$GOs,",")
gene2go <- data.frame(GID = rep(gterms$query,
                                   times = sapply(all_go_list, length)),
                         GO = unlist(all_go_list),
                         EVIDENCE = "IEA")
  						 						 
1.5 #提取基因与KEGG的对应
gene2ko <- egg %>%
  dplyr::select(GID = query, KO = KEGG_ko) %>%
  na.omit()	
  
gene2ko$KO[which(gene2ko$KO=="-")]=NA 
gene2ko <- na.omit(gene2ko) 

1.6 #构建了一个函数，利用KEGG的json来获取ko与通路以及K与ko的对应关系
if(!file.exists('kegg_info.RData')){

  library(jsonlite)
  library(purrr)
  library(RCurl)

  update_kegg <- function(json = "ko00001.json",file=NULL) {
    pathway2name <- tibble(Pathway = character(), Name = character())
    ko2pathway <- tibble(Ko = character(), Pathway = character())

    kegg <- fromJSON(json)

    for (a in seq_along(kegg[["children"]][["children"]])) {
      A <- kegg[["children"]][["name"]][[a]]

      for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
        B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 

        for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
          pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]

          pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
          pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
          pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))

          kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]

          kos <- str_match(kos_info, "K[0-9]*")[,1]

          ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
        }
      }
    }

    save(pathway2name, ko2pathway, file = file)
  }

  update_kegg(json = "ko00001.json",file="kegg_info.RData")

}
 
load(file = "kegg_info.RData")

1.7 #利用`gene2ko`与`ko2pathway`将基因与pathway对应起来
# 在运行数据框合并前，需要做到两个数据框的列名是对应的，并且将原来gene2ko中的ko修改一下
colnames(ko2pathway)=c("KO",'Pathway')
library(stringr)
gene2ko$KO=str_replace(gene2ko$KO,"ko:","")

# 合并代码是：
gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "KO") %>% 
  dplyr::select(GID, Pathway) %>%
  na.omit()
  
#根据NCBI Taxonomy 查询物种的Taxonomy，例如要查sesame
tax_id = "3750"
genus = "Malus" 
species = "domestica"

library(AnnotationForge)
makeOrgPackage(gene_info=gene_info,
               go=gene2go,
               ko=gene2ko,
               maintainer='yjlin <yjlin@psc.ac.cn>',
               author='yjlin',
               pathway=gene2pathway,
               version="0.0.1",
               outputDir = ".",
               tax_id="3750",
               genus="Malus",
               species="domestica",
               goTable="go")
			   
1.9 #像安装R包一样，来安装物种注释包
install.packages("./org.Mdomestica.eg.db", repos=NULL, type="source")  #, type="source"
# 然后加载
library(org.Mdomestica.eg.db) 
library(clusterProfiler)

pathway2gene <- AnnotationDbi::select(org.Mdomestica.eg.db, 
                                      keys = keys(org.Mdomestica.eg.db), 
                                      columns = c("Pathway","KO")) %>%
  na.omit() %>%
  dplyr::select(Pathway, GID)
  

2. #进行KEGG富集分析
load("kegg_info.RData")

"""
x = enricher(glist, TERM2GENE=pathway2gene, TERM2NAME=pathway2name)
head(summary(x),1)
"""

#input file
setwd("J:/C盘/Project/psc-langlab/task/apple/wangxiaofei/2021年8月19日-DREB2A/ko/kegg")
glist <- read.table("DREB2A.IDR_All.geneid.txt",header=FALSE) #target gene list/DEGs
glist <- glist[,1]
geneList <- read.table("DREB2A.IDR_All.geneid.txt",header=FALSE) #background gene list
geneList <- geneList[,1]

#KEGG pathway 富集
ekp <- enricher(glist, universe=geneList,
                TERM2GENE = pathway2gene, 
                TERM2NAME = pathway2name, 
                pvalueCutoff = 1, 
                qvalueCutoff = 1,
                pAdjustMethod = "BH",
                minGSSize = 1)

ekp_results <- as.data.frame(ekp)

barplot(ekp, showCategory=20,color="pvalue",
        font.size=10)
dotplot(ekp)

emapplot(ekp)
