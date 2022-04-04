
library(data.table)
library(ggplot2)

l <- 80

dt2 <- fread("../../Data/true_module_term_bp.csv")
# 62 out of 80 modules are enriched in biological processes 2021
enrmodbp <- table(dt2$terms != 0)[[2]]

dt3 <- fread("../../Data/shuffle_enr_term_bp.csv")

dtg2 <- data.table(num_enr_terms = colSums(dt3))
dtg2[, avg := num_enr_terms / l]

g2 <- ggplot(dtg2, aes(x = avg)) +
    geom_histogram(aes(y = ..density..), alpha = 0.2,
        position = "identity",
        #fill = NA,
        bins = 30) +
    geom_density(aes(avg),
        color = "cyan3", fill = "cyan1",
        alpha = 0.1,
        size = 1) +
    geom_vline(xintercept = sum(dt2) / l, color = "blue", linetype = "dashed",
        size = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
    scale_x_continuous(breaks = c(seq(15, 35, 5), round(sum(dt2) / l))) +
    xlab("Average number of enriched biological processes per gene module") +
    ylab("Density")

dt <- fread("../../Data/enr_num_bp.csv")

g1 <- ggplot(dt, aes(x = enr100)) +
    geom_histogram(aes(y = ..density..), alpha = 0.2,
        position = "identity",
        #fill = NA,
        binwidth = 1) +
    geom_density(aes(enr100),
        color = "cyan3", fill = "cyan1",
        alpha = 0.1,
        size = 1) +
    geom_vline(xintercept = enrmodbp, color = "blue", linetype = "dashed",
        size = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
    scale_x_continuous(breaks = c(45, 50, 55, 60, 62),
        name = expression("Number of gene modules with ">=" 1 enriched biological processes")) +
    ylab("Density")

pdf("../../results/FigureS9.pdf", width = 10)
print(g1)
print(g2)
dev.off()







