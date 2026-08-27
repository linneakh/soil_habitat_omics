### FUNCTIONS FOR PCA and nmds

list_of_shapes <- c(15,17,18,3,7,8,9,10,4,11,12,13,14)

col_list2 = colors = c("black","gray","coral4","coral",
                       "chartreuse3", "darkseagreen1","blue","lightblue",
                       "yellow4", "yellow","darkorchid", "plum2",
                       "darkred", "darksalmon","green4",
                       "greenyellow","orange",
                       "moccasin",   "hotpink4")
col_list_zone = c("green", "blue", "orange",  "red")
col_list_treatment = c("darkred", "darkgreen")



make_nmds_plot <- function(nmds_object, Group1, Group2, col_list=NULL){
  grp1 <- enquo(Group1)
  grp2 <- enquo(Group2)
  nmds_plot <- ggplot(data=nmds_object, mapping = aes(x=NMDS1, y=NMDS2, col=!!grp1, shape=!!grp2)) +
    geom_point(size=size/3, show.legend = TRUE) +
    theme_linedraw(base_size = size) + labs(x= "NMDS1", y="NMDS2") +
    scale_color_manual(values = c(col_list)) +
    scale_shape_manual(values= list_of_shapes) +
    theme( legend.text = element_text(size=size+3, face="bold"),
           legend.title = element_blank(),
           legend.key.size = unit(0.6, "cm"),
           legend.key.width = unit(0.6,"cm"),
           #legend.position = "bottom",
           axis.title.x = element_text(size=size+3,face="bold"),
           axis.title.y = element_text(size=size+3,face="bold"),
           plot.title = element_text(size=size+3,face="bold"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank()) 
  return(nmds_plot)
}



