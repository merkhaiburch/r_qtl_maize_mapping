# ----------------------------------
# Merritt Burch
# mbb262@cornell.edu
# 2019-07-24
# Script to plot Ga1-s results and 
# Overlay the two plots
# ----------------------------------

# Must run B73xM162 and B73xKy21 files first to get files

# rbind the two dataframes together
out.ehk_ky21$id <- rep("Ky21", nrow(out.ehk_ky21))
out.ehk_m162w$id <- rep("M162w", nrow(out.ehk_m162w))
results <- rbind(out.ehk_ky21, out.ehk_m162w)

# Send results to result folder
setwd("C:/Users/merri/git_projects/ga1s/results")

# Plot alone
# pdf("all.pdf", width = 7.5, height = 3.5)
# ggplot(data=out.ehk_ky21, aes(x=pos, y = lod)) + 
#   geom_line(size = 0.9, alpha = .8) + 
#   facet_wrap(~chr, nrow = 1, scales = "free_x") +
#   theme_classic() +
#   theme(strip.text = element_text(face = "bold", size = 15),
#         axis.text.y = element_text(size = 10, colour = "black"),
#         axis.text.x = element_text(angle= 60, hjust=0.8, vjust = 1, size = 10, colour = "black"),
#         strip.background = element_rect(colour = "white"), 
#         axis.line = element_line(color = "black", size = 1),
#         panel.spacing = unit(.2, "lines"),
#         axis.ticks = element_line(colour = "black")) +
#   geom_hline(yintercept = 3, color = "red") +
#   scale_y_continuous(limits = c(0, 7.5), expand = c(0.02, 0)) +
#     labs(x = "Chromosome Position",
#          y = "LOD Score",
#          color = "",
#          linetype = "")
# dev.off()

# Plot together
pdf("both_Ky21_and_M162w.pdf", width = 6.5, height = 4)
ggplot(data = results, aes(x=pos, y = lod, group = id)) + 
  geom_line(size = 0.8, alpha = 0.9, aes(color = id)) + 
  facet_wrap(~chr, nrow = 1, scales = "free_x") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 15),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle= 60, hjust=0.8, vjust = 1, size = 10, colour = "black"),
        strip.background = element_rect(colour = "white"), 
        axis.line = element_line(color = "black", size = 1),
        panel.spacing = unit(.2, "lines"),
        axis.ticks = element_line(colour = "black"),
        legend.position= "top") +
  geom_hline(yintercept = 3, color = "red") +
  scale_color_manual(values=c("#6b6b6b", "#000000")) +
  # scale_color_grey() +
  scale_y_continuous(limits = c(0, 7.5), expand = c(0.02, 0)) +
  labs(x = "Chromosome Position",
       y = "LOD Score",
       color = "",
       linetype = "")
dev.off()



# ---------------------------
#   Plot additivity effects
# ---------------------------

# Generate the data
tmp_ky21 <- sim.geno(ga_ky21, map.function = "haldane")
effectscans_ky21 <- effectscan(tmp_ky21, chr=1:10, draw = F)

tmp_m162 <- sim.geno(ga, map.function = "haldane")
effectscans_m162 <- effectscan(tmp_m162, chr=1:10, draw = F)

# Add identifiers then Combine the data
effectscans_ky21$id <- rep("Ky21", nrow(effectscans_ky21))
effectscans_m162$id <- rep("M162w", nrow(effectscans_m162))
additive_effect <- rbind(effectscans_ky21, effectscans_m162)

# Whole genome (single plot)
# ggplot(data=effectscans, aes(x=pos, y = a)) + 
#   geom_line(size = 0.8, alpha = .8) + 
#   geom_hline(yintercept = 0, linetype = "solid") +
#   facet_wrap(~chr, nrow = 1, scales = "free_x") +
#   theme_classic() +
#   theme(strip.text = element_text(face = "bold", size = 15),
#         axis.text.y = element_text(size = 15, colour = "black"),
#         axis.text.x = element_text(angle=45,hjust=1, size = 12, colour = "black"),
#         strip.background = element_rect(colour = "white"), 
#         axis.line = element_line(color = "black", size = 1),
#         panel.spacing = unit(.2, "lines")) +
#   scale_y_continuous(limits = c(-.35, .35), expand = c(0, 0)) +
#   labs(x = "Chromosome Position",
#        y = "Additive Effect",
#        color = "",
#        linetype = "")


# Whole genome (both RILs)
pdf("both_Ky21_and_M162w_effects.pdf", width = 6.8, height = 3)
ggplot(data=additive_effect, aes(x=pos, y = a, group = id)) + 
  geom_line(size = 0.8, alpha = .8, aes(color = id)) + 
  geom_hline(yintercept = 0, linetype = "solid") +
  facet_wrap(~chr, nrow = 1, scales = "free_x") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 15),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle=60,hjust=0.8, vjust = 1, size = 10, colour = "black"), 
        axis.line = element_line(color = "black", size = 1),
        panel.spacing = unit(.2, "lines"),
        axis.ticks = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")  +
  scale_color_manual(values=c("#6b6b6b", "#000000")) +
  scale_y_continuous(limits = c(-.55, .55), expand = c(0, 0)) +
  labs(x = "",
       y = "Additive Effect",
       color = "",
       linetype = "")
dev.off()
