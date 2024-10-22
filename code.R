library(ggplot2)
library(limma)
library(annotables)
library(Biobase)
library(tidyverse)

g <- glimpse  
h <- head        
s <- summary       


ps <- read_tsv("GSE47944_ps_tt_gel_raw.txt")


replace_letter <- function(vec) {
        for (i in 1:length(vec)) {
                string <- vec[i]
                string_vec <- unlist(strsplit(string, ""))
                idx <- which("|" == string_vec)
                string_vec[idx] <- " "
                vec[i] <- paste(string_vec, collapse = "")
                
        }
        return(vec)
}



ps1 <- ps %>%
        mutate(Gene = replace_letter(Gene)) %>%
        separate(Gene, c("ID", "symbol"), sep = " ")

# ps_f: feature matrix
ps_2 <- ps1 %>%
        transmute(ensgene = ID, 
                  symbol = symbol) %>%
        left_join(grch38, by = "symbol", suffix = c("_ps", "_db")) %>%
        filter(!is.na(entrez))

ps_3 <- ps_2[, c("symbol", "entrez", "chr")] %>%
        filter(!duplicated(symbol)) %>%
        mutate(entrez = as.character(entrez))

ps_4 <- as.matrix(ps_3) 
rownames(ps_4) <- ps_3$symbol

ps_f <- as.data.frame(ps_4)



# ps_x: expression matrix 


ps_5 <-  subset(ps1, symbol %in% ps_f$symbol) %>%
        filter(!duplicated(symbol))

ps_6 <- ps_5[, colnames(ps)[-1]] %>%
        as.matrix()

rownames(ps_6) <- ps_f$symbol
ps_x <- ps_6


# ps_p: phenotype matrix

trt <- rep(c("Untreated", "DMSO", "Agonist", "Antagonist"), times = 21)

ps_9 <- data.frame(sample = colnames(ps_x)) %>%
        mutate(disease = ifelse(str_detect(sample, "K"), 
                                "Psoriasis", 
                                "Normal"))

ps_9$treatment <- trt


colnames(ps_x) <- ps_9$sample
ps_10 <- as.matrix(ps_9)
rownames(ps_10) <- ps_9$sample 
ps_p <- as.data.frame(ps_10)[, -1]


############################## Create ExpressionSet ################################

# eset: ExpressionSet object
eset <- ExpressionSet(assayData = ps_x,
                      phenoData = AnnotatedDataFrame(ps_p),
                      featureData = AnnotatedDataFrame(ps_f))




################################ pre-precessing ##################################

# Log transform
exprs(eset) <- log(exprs(eset))
plotDensities(eset, 
              group = pData(eset)[, "disease"],
              legend = "topright", 
              main = "Log-Transformed Count Distribution without Normalization")


# normalization
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
plotDensities(eset, 
              group = pData(eset)[, "disease"],
              legend = "topright", 
              main = "Log-Transformed Count Distribution with Normalization")
abline(v = 6, col = "blue")


# filtering 
keep <- which(rowMeans(exprs(eset)) > 6)
eset <- eset[keep, ]
plotDensities(eset, 
              grou = pData(eset)[, "disease"],
              legend = "topright",
              main = "Count Distribution after Filtering")






################################ Sample Inspection ################################

# hierarchical clustering & heatmap
cor_mat <- cor(exprs(eset))
library(pheatmap)
pheatmap(cor_mat,
         annotation = pData(eset)[, c("disease", "treatment")],
         main = "Hierarchical Clustering")

# PCA 
plotMDS(eset, 
        labels = pData(eset)[, "disease"],
        gene.selection = "common",
        main = "PCA by Disease Status")


plotMDS(eset, 
        labels = pData(eset)[, "treatment"],
        gene.selection = "common",
        main = "PCA by Treatment")


################################### Modeling ###################################

group <- with(pData(eset), paste(disease, treatment, sep = "_"))
group <- factor(group)

design <- model.matrix(~ 0 + group)

colnames(design) <- levels(group)

biological_replicates <- colSums(design) 

# cm: contrast matrix

cm <- makeContrasts(disease_untreated = Psoriasis_Untreated - Normal_Untreated,
                    disease_dmso = Psoriasis_DMSO - Normal_DMSO,
                    disease_agonist = Psoriasis_Agonist - Normal_Agonist,
                    disease_antagonist = Psoriasis_Antagonist - Normal_Antagonist,
                    normal_dmso = Normal_DMSO - Normal_Untreated,
                    normal_agonist = Normal_Agonist - Normal_DMSO,
                    normal_antagonist = Normal_Antagonist - Normal_DMSO,
                    psoriasis_dmso = Psoriasis_DMSO - Psoriasis_Untreated,
                    psoriasis_agonist = Psoriasis_Agonist - Psoriasis_DMSO,
                    psoriasis_antagonist = Psoriasis_Agonist - Psoriasis_DMSO,
                    interaction = (Normal_DMSO - Normal_Untreated) - 
                            (Psoriasis_DMSO - Psoriasis_Untreated),
                    levels = design)


fit <- lmFit(eset, design)
fit1 <- contrasts.fit(fit, contrasts = cm)
fit2 <- eBayes(fit1)
res <- decideTests(fit2)

# additional data cleaning for plotting 
res_summary0 <- as.data.frame(summary(res)) %>%
        rownames_to_column()

res_summary <- res_summary0[, -1] %>%
        spread(Var1, Freq)

names(res_summary) <- c("Group", "Down-regulated", "Insignificant", "Up-regulated")
res_summary1 <- gather(res_summary, Gene_Category, Number_of_Genes, -Group) %>%
        mutate(Gene_Category = factor(Gene_Category, 
                                      levels = c("Down-regulated", 
                                                 "Insignificant", 
                                                 "Up-regulated")))

res_summary1$Number_of_Genes <- as.numeric(res_summary1$Number_of_Genes)



############################### Result Inspection ###############################

get_results_table <- function(cf) {
        topTable(fit2, 
                 coef = cf,
                 number = nrow(fit2), 
                 sort.by = "none")
}

stats_disease_untreated <- get_results_table("disease_untreated")

stats_disease_dmso <- get_results_table("disease_dmso")

stats_disease_agonist <- get_results_table("disease_agonist")

stats_disease_antagonist <- get_results_table("disease_antagonist")

stats_normal_dmso <- get_results_table("normal_dmso")

stats_normal_agonist <- get_results_table("normal_agonist")

stats_normal_antagonist <- get_results_table("normal_antagonist")

stats_psoriasis_dmso <- get_results_table("psoriasis_dmso")

stats_psoriasis_agonist <- get_results_table("psoriasis_agonist")

stats_psoriasis_antagonist <- get_results_table("psoriasis_antagonist")

stats_interaction <- get_results_table("interaction")

# P value data frame 
pval_df <- data.frame(disease_untreated = stats_disease_untreated$P.Value,
                      disease_dmso = stats_disease_dmso$P.Value,
                      disease_agonist = stats_disease_agonist$P.Value,
                      disease_antagonist = stats_disease_antagonist$P.Value,
                      normal_dmso = stats_normal_dmso$P.Value,
                      normal_agonist = stats_normal_agonist$P.Value,
                      normal_antagonist = stats_normal_antagonist$P.Value,
                      psoriasis_dmso = stats_psoriasis_dmso$P.Value,
                      psoriasis_agonist = stats_psoriasis_agonist$P.Value,
                      psoriasis_antagonist = stats_psoriasis_antagonist$P.Value,
                      interaction = stats_interaction$P.Value) %>%
        gather(Contrast, 
               P.Value) %>%
        mutate(Contrast = factor(Contrast,
                                 levels = c("disease_untreated",
                                            "disease_dmso",
                                            "disease_agonist",
                                            "disease_antagonist",
                                            "normal_dmso",
                                            "normal_agonist",
                                            "normal_antagonist",
                                            "psoriasis_dmso",
                                            "psoriasis_agonist",
                                            "psoriasis_antagonist",
                                            "interaction")))


################################## Data Cleaning for Presentation ############################
gene_symbols <- fit2$genes[, "symbol"]

clean_data <- function(df, label_threshold) {
        df %>% mutate(LogOdds = -log10(adj.P.Val),
                      Category = case_when(adj.P.Val < 0.05 & logFC > 0 ~ "Up-regulated",
                                           adj.P.Val < 0.05 & logFC < 0 ~ "Down-regulated",
                                           adj.P.Val > 0.05 ~ "Insignificant"),
                      Label = ifelse(LogOdds > label_threshold, symbol, "")) %>%
                mutate(Category = factor(Category, 
                                         levels = c("Up-regulated",
                                                    "Insignificant",
                                                    "Down-regulated")))
        
}

# data frames for volcano plots
stats_psoriasis_dmso_1 <- clean_data(stats_psoriasis_dmso, 25)
stats_normal_dmso_1 <- clean_data(stats_normal_dmso, 12.5)
stats_disease_untreated_1 <- clean_data(stats_disease_untreated, 4)

# data frames for a heat map
stats_disease_untreated_2 <- subset(stats_disease_untreated, 
                                    adj.P.Val < 0.05)
untreated_samples <- subset(pData(eset), treatment == "Untreated")


count_disease_untreated <- exprs(eset)[stats_disease_untreated_2$symbol, ]
count_disease_untreated1 <- count_disease_untreated[, rownames(untreated_samples)]

# data frame for affected gene list
disease_fc <- subset(stats_disease_untreated, 
                     adj.P.Val < 0.05) %>%
        arrange(desc(logFC)) 
disease_fc1 <- rbind(disease_fc[1:10, ], 
                     disease_fc[451:460, ]) %>%
        mutate(Category = ifelse(logFC > 0, 
                                 "Up-regulated", 
                                 "Down-regulated"),
               Category = factor(Category,
                                 levels = c("Up-regulated", 
                                            "Down-regulated")))


###################################### GSEA ####################################


# data cleaning for GO GSEA
disease_untreated_gene_list <- disease_fc$logFC
names(disease_untreated_gene_list) <- disease_fc$symbol
disease_untreated_gene_list1 = sort(disease_untreated_gene_list, decreasing = TRUE)

library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

# creating a GO object
disease_untreated_GO <-
        gseGO(geneList = disease_untreated_gene_list1, 
              ont = "ALL", 
              keyType = "SYMBOL", 
              nPerm = 10000, 
              minGSSize = 3, 
              maxGSSize = 800, 
              pvalueCutoff = 0.05, 
              verbose = TRUE, 
              OrgDb = org.Hs.eg.db, 
              pAdjustMethod = "none")

# data cleaning for KEGG GSEA
disease_untreated_gene_list2 <- disease_fc$logFC
names(disease_untreated_gene_list2) <- disease_fc$entrez
disease_untreated_gene_list3 <- sort(disease_untreated_gene_list2, decreasing = TRUE)

# creating a KEGG object
disease_untreated_KEGG <-
        gseKEGG(geneList     = disease_untreated_gene_list3,
                organism     = "hsa",
                nPerm        = 10000,
                minGSSize    = 3,
                maxGSSize    = 800,
                pvalueCutoff = 0.05,
                pAdjustMethod = "none",
                keyType       = "ncbi-geneid")







#################################### Plotting #####################################

contrasts_summary_plot <- 
        ggplot(res_summary1, 
               aes(x = Number_of_Genes,
                   y = Group,
                   color = Gene_Category)) + 
        geom_point(alpha = 0.5,
                   size = 3) + 
        theme_bw() +
        ylab("Contrasts") + 
        xlab("Affected Gene Number") + 
        ggtitle("Summary of Gene Expression Change")

pval_plot <- 
        ggplot(pval_df,
               aes(x = P.Value,
                   fill = Contrast,
                   color = Contrast)) + 
        geom_density(alpha = 0.3) +
        theme_bw() + 
        xlab("P-Value") + 
        ylab("Proportion") + 
        ggtitle("Distribution of P-Values over Contrasts")


normal_psoriasis_heatmap <- 
        pheatmap(count_disease_untreated1,
                 annotation = untreated_samples[, c("disease", "treatment")],
                 main = "Normalized Gene Expression Profile in Normal vs Psoriasis Skin")

gene_rank_plot <- 
        ggplot(disease_fc1, 
               aes(x = logFC,
                   y = reorder(symbol, logFC),
                   fill = Category)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values=c("red", "blue")) +
        theme_bw() +
        xlab("Log2 Fold Change") +
        ylab("Genes") + 
        ggtitle("10 Most Up- and Down-regulated Genes") +
        scale_x_continuous(limits = c(-5, 5))


library(ggrepel)

volcano_plot <- function(df, tit) {
        ggplot(df,
               aes(x = logFC,
                   y = LogOdds,
                   color = Category,
                   label = Label)) +
                geom_point(alpha = 0.3) + 
                theme_bw() + 
                scale_x_continuous(limits = c(-6, 6)) + 
                xlab("Log2 Fold Change") +
                ylab("-Log10 Adjusted P-value") + 
                ggtitle(tit) +
                scale_color_manual(values = c("red", "#999999", "blue")) + 
                geom_text_repel(color = "black")
}

volcano_psoriasis_dmso <- volcano_plot(stats_psoriasis_dmso_1,
                                       "Gene Expression in Psoriasis Skin under DMSO vs Untreated Condition")

volcano_normal_dmso <- volcano_plot(stats_normal_dmso_1,
                                    "Gene Expression in Normal Skin under DMSO vs Untreated Condition")

volcano_disease_untreated <- volcano_plot(stats_disease_untreated_1,
                                          "Gene Expression in Psoriasis vs Normal Skin under Untreated Condition")



library(gridExtra)

grid.arrange(volcano_normal_dmso,
             volcano_psoriasis_dmso,
             volcano_disease_untreated,
             nrow = 3)

# GSEA plots



disease_untreated_GO_dot <- 
        dotplot(disease_untreated_GO, 
                showCategory = 10, 
                split=".sign") + 
        facet_grid(.~.sign) + 
        ggtitle("Biological Events in Psoriasis vs Normal Skin: GO") + 
        theme(axis.text.y = element_text(size = 12))


disease_untreated_GO_ridge <-
        ridgeplot(disease_untreated_GO) + 
        labs(x = "Enrichment Distribution") + 
        ggtitle("Fold Change Distribution") + 
        theme(axis.text.y = element_text(size = 11)) 

disease_untreated_GO_emap <- 
        emapplot(disease_untreated_GO, 
                 showCategory = 10)


disease_untreated_GO_cnet <- 
        cnetplot(disease_untreated_GO, 
                 categorySize = "enrichmentScore", 
                 foldChange = disease_untreated_gene_list1, 
                 showCategory = 3)

disease_untreated_KEGG_dot <- 
        dotplot(disease_untreated_KEGG, 
                showCategory = 10, 
                title = "Biological Events in Psoriasis vs Normal Skin: KEGG" , 
                split=".sign") + 
        facet_grid(.~.sign) +
        theme(axis.text.y = element_text(size = 12))

disease_untreated_KEGG_ridge <-
        ridgeplot(disease_untreated_KEGG) + 
        labs(x = "Enrichment Distribution") + 
        ggtitle("Fold Change Distribution") + 
        theme(axis.text.y = element_text(size = 11)) 
