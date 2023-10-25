



phenotypePlot <-
  function(d, max.y,max.x, suggestive.line, significant.line, 
           size.x.labels=9, size.y.labels=9,  sort.by.value=F, sort.by.category.value=F,
           #annotation
           annotate.angle=0, annotate.size=5, 
           annotate.list, 
           #labels
           x.group.labels=TRUE, 
           sizes=F, direction=F, point.size=3.5,
           #plot characteristics
           use.color=T,
           color.palette,
           title= paste0("Phenotype Plot ", date()),
           x.axis.label="Phenotypes",
           y.axis.label="Values",
           y.axis.interval=5) {
	
    #Check to ensure the input contains columns phenotype and value
    if ( sum(c("phenotype","value") %in% names(d))<2 ) stop("Data input must contain columns phenotype and value.")
    
    #Check for conflicting color commands
    if(!use.color & !missing(color.palette)) stop("You requested no color, but provided a color palette.")
    #Check for color
    if(use.color) {
      if(missing(color.palette)) {
        if(!("color" %in% names(d))) stop("You requested color, but did not provide a color attribute in d or color.palette")
        else if(class(d$color)=="factor") warning("The color attribute is a factor and no color palette is provided: R default color scheme will be used. Convert color to character if it contains color names or codes")
      }
      #Set up color if using a palette
      if(!missing(color.palette)&!length(d$color)) {
        if(length(d$groupnum)) d$color=d$groupnum
        else stop("You requested use.color, but d$color or d$groupnum were not provided to distinguish groups.")
      }
    } else {
      #Set the colors to all black if no color is requested
      d$color="#000000"
    }
    
    #Check for point sizing/direction information if requested
    if(sizes&!length(d$size)) stop("You requested size information, but did not provide d$size")
    if(direction&!length(d$direction)) stop("You requested direction information, but did not provide d$direction")
    #Set point size if sizes are not requested
    if(!sizes) d$size=point.size
    
    #One cannot sort by value and have x.group.labels
    if (sort.by.value & x.group.labels) stop("One cannot sort universally by value and have x group labels. Try sort.by.category.value or not labeling the groups.")
    #Check for group information if requested
    if((sort.by.category.value|x.group.labels)&!("groupnum" %in% names(d))) {
      stop("Requested group information, but did not provide d$groupnum")
    }
    
    #If no group info, just assign all to group 0
    if(!("groupnum"%in%names(d))) { 
      d$groupnum=0
    }
    
    #Remove items with NA values or phenotypes
    d=d[!(is.na(d$phenotype)|is.na(d$value)),]
    
    #Sort by the phenotype
    d=d[order(d$phenotype),]
        
    #Set the maximum x value to fit all phenotypes if not specified
    if(missing(max.x)) max.x = length(unique(d$phenotype))
    
    
    #Create the list of phenotypes, finding the best values for each phenotype
    phenotypes=aggregate(value ~ phenotype + groupnum, d,FUN=max)
    #Remove the least significant phenotypes; only has an effect if max.x was specified.
    phenotypes=phenotypes[order(phenotypes$value, decreasing=T),][1:min(nrow(phenotypes),max.x),]
    
    #If the user requested sorting by values
    if (sort.by.value) {
      #Sort by values
      phenotypes=phenotypes[order(phenotypes$value, decreasing=T),]
    }  else if (sort.by.category.value) {
      #If the user requested sorting by values within each group.
      #Sort by group and then value
      phenotypes=phenotypes[order(-phenotypes$groupnum,phenotypes$value, decreasing=T),]
    } else {
      #Restore the phenotype sorting order if no other sorting was requested
      #TODO: What about non-numeric phenotypes? Need this to sort 008 etc.
      phenotypes=phenotypes[order(phenotypes$groupnum,phenotypes$phenotype),]
    }
    
    phenotypes$seq = 1:nrow(phenotypes)      

    
    #Limit to phenotype and seq, as they are the only relevant columns
    #Include value as min.value for annotation purposes
    phenotypes=phenotypes[,c("phenotype","seq","value")]
    names(phenotypes)[3]="min.value"
    
    #Add sequence information
    d=inner_join(phenotypes,d,by="phenotype")
    d=d[order(d$seq),]
    
    #Define the max y axis value if not provided
    if(missing(max.y)) max.y=ceiling(max(d$value))
    
    if (x.group.labels) {
      labels= summarize(group_by(d, groupnum), tick=mean(unique(seq)),label=as.character(groupnum[1]))
      labels=labels[order(labels$tick),]
    }
    
    if(missing(color.palette)) {
      color.palette = unique(d[order(d$seq),]$color)
      names(color.palette)=color.palette
    } else {
      names(color.palette) = unique(d[order(d$seq),]$color)
    }
    

      #Generate the inital plot
	  
      plot=ggplot(d,ylab=y.axis.label,xlab=x.axis.label)
	  
      #Include lines for significance thresholds
      if (!missing(suggestive.line)&!is.na(suggestive.line)) plot=plot+geom_hline(yintercept=suggestive.line,colour="black", alpha=I(1/3),linewidth=1, linetype="dashed")
      if (!missing(significant.line)&!is.na(significant.line)) plot=plot+geom_hline(yintercept=significant.line,colour="#00AE00",alpha=I(1/3),linewidth=1)
    #   plot = plot+geom_hline(yintercept=7.301, colour="#191970",alpha=I(1/3),size=1) ## Ajout d'une ligne au genome-wide threshold
      
      plot=plot+aes(seq,value,size=size,colour=color)
      if(!sizes) plot=plot+scale_size(range=c(point.size,point.size),guide="none")
      #Add points
      plot=plot+geom_point(alpha =0.5) 
      
      #Color as defined 
      plot = plot + scale_colour_manual(values= color.palette, guide="none") 
      
      #If label the X axis with the groups if requested
      if (x.group.labels) {
      
        plot=plot+scale_x_continuous(name=x.axis.label, limits=c(1,max.x), breaks=labels$tick, labels=labels$label, expand=c(.01,0))
        
      } else {
        plot=plot+scale_x_continuous(name=x.axis.label, limits=c(1,max.x), breaks=c(-100), labels=c(""), expand=c(.015,0))
      }
      
      #Set the Y scale and labels
      #plot=plot+scale_y_continuous(y.axis.label, limits=c(0,max.y), breaks=seq(0,max.y,1), expand=c(0,0.02))
	  plot=plot+scale_y_continuous(y.axis.label, limits=c(0,max.y), breaks=seq(0,max.y,y.axis.interval), expand=c(0,0.02))
      
      #Set the default theme
      plot=plot+theme(
	    
		panel.background = element_rect(fill = "white"),
        #panel.background=element_blank(),
		plot.margin = margin(0.5,1.5,0.5,0.5,"cm"),
		plot.background = element_rect(fill = "white",colour = "white",linewidth = 1),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(size=size.x.labels, colour="black", angle=-50, hjust=0, vjust=1, face="bold"), 
        axis.text.y=element_text(size=size.y.labels, colour="black", face = "bold"), 
        axis.line =element_line(colour="black"),
        axis.ticks=element_line(colour="black")
      ) 
      
 
    #Hide the legend by default
    plot = plot+theme(legend.position = "none")
    
    #Add OR information
    if(sizes){			
      plot= suppressWarnings(plot + theme(legend.position="right") +
        scale_size("Size", range = c(point.size, 2*point.size)))
    }
    if(direction){
      plot=plot+aes(shape = factor(direction), fill=color) + 
        scale_shape("Direction", solid = TRUE) + scale_shape_manual(values=c(25,24)) + 
        scale_fill_manual(values= color.palette, guide="none") 
    }		
    #If annotation present, start the definitions, otherwise skip it
    plot = plot + ggrepel::geom_text_repel(aes(label=description),colour="black",data=d[d$annotate,],size=annotate.size,angle=annotate.angle, direction = "both")   # Annotation originale
    
    #Add the title
    plot=plot+labs(title=title) + theme(title=element_text(size=15))
    
  }


#' This function take as input a data.frame of primary MR analysis and return a balloon plot. credit Nicholas.
#'
#' @param dat_plot a data.frame containing the column "exposure" containing all the name of your exposure (or outcomes) as factor. the names will appear in the order you specified in the leevels of your factor from the bottom to up.
#' containing the column "outcome"  containing all the name of your outcomes (or exposures) as factor. the names will appear in the order you specified in the leevels of your factorfrom left to right.
#' containing the column "log10_pval" a numeric column containing the log(p-value) or any measure that will be scaled with the size of the circle 
#' The higher the measure the larger the circle. Containing the column "shape_point" which can take two values "circle" or "cross" which will become cross.
#' containing the column "z_score" a numeric column containing the z-score or beta which will be attributed to a colour gradient.
#' containing the column "Category"  a factor column the name of the category. the names will appear in the order you specified in the levels of your factor from left to right on top of the graph
#' If the function throws an error it is because you did not specified dat_plot accordingly so read this section carefully.
#' @param subdivide_by_category if TRUE will divide by category. If FALSE won't
#' @param name_pval the name of the legend that is mapped to log10_pval
#' @param name_zscore the name of the legend that is mapped to z_score
#' @return A Ballon plot of the primary MR analysis.
#' @export
plot_balloon2 <- function(dat_plot, subdivide_by_category = TRUE, 
                          name_pval = expression(-Log[10](P)), name_zcore = "Z-score") {
  stopifnot(c("exposure", "outcome", "log10_pval", "shape_point", "z_score", "Category") %in% colnames(dat_plot))
  dat_plot_rond = dat_plot
  dat_plot_rond[shape_point == "cross", z_score := NA]
  dat_plot_rond[shape_point == "cross", log10_pval := NA]
  # dat_plot_rond$shape_point = sapply(dat_plot_rond$pval, FUN = function(x) {ifelse( x < bonferroni_threshold, "rond", NA)})
  
  
  dat_plot_croix = dat_plot_rond
  dat_plot_croix$shape_point <- ifelse(is.na(dat_plot$z_score), "Non-significant", NA)
  
  # dat_plot_croix$shape_point = sapply(dat_plot_croix$pval, FUN = function(x) {ifelse( x >= bonferroni_threshold, "Non-significant", NA)})
  # 
  ordered_category = dat_plot %>% dplyr::arrange(Category) %>% dplyr::select(Category) %>% dplyr::distinct()
  breaks_category = sapply(unique(dat_plot$Category), FUN = function(x){dim(dat_plot %>% dplyr::select(outcome , Category) %>% dplyr::distinct() %>% dplyr::filter(Category ==  !!x))[1]})
  
  for (i in 1:length(ordered_category$Category[!is.na(ordered_category$Category)])) {
    if (i == 1)
    {
      list_breaks = as.numeric(breaks_category[i])
      label_breaks = as.numeric(breaks_category[i])/2
    } else {
      list_breaks = c(list_breaks, (data.table::last(list_breaks) + as.numeric(breaks_category[i])))
      label_breaks = c(label_breaks, (list_breaks[i - 1] + (as.numeric(breaks_category[i])/2 + 0.5)))
    }
  }
  
  balloon_plot = ggplot2::ggplot() +
    ggplot2::geom_point(data = dat_plot_croix, ggplot2::aes(x = outcome, y = exposure, shape = factor(shape_point)), size = 2, color = "gray20")
  
  when_annotate = seq(from = 1 , to  = length(ordered_category$Category[!is.na(ordered_category$Category)]), by = 2)
  
  cpt = 0
  for (i in c(1:length(ordered_category$Category[!is.na(ordered_category$Category)]))) {
    cpt = cpt + 1
    if (i == 1)
    {
      list_breaks_background = as.numeric(breaks_category[i]) + 0.5
      balloon_plot = balloon_plot +
        ggplot2::annotate(
          "rect",
          xmin = 0,
          xmax = list_breaks_background[i],
          ymin = 0,
          ymax = length(unique(dat_plot$exposure)) + 0.5,
          fill = "gray5",
          alpha = 0.1
        )
    } else {
      list_breaks_background = c(list_breaks_background, (data.table::last(list_breaks_background) + as.numeric(breaks_category[i])))
      if (cpt %in% when_annotate)
      {
        balloon_plot = balloon_plot +
          annotate(
            "rect",
            xmin = list_breaks_background[i - 1],
            xmax = list_breaks_background[i],
            ymin = 0,
            ymax = length(unique(dat_plot$exposure)) + 0.5,
            fill = "gray5",
            alpha = 0.1
          )
      }
    }
  }
  
  
  balloon_plot <- balloon_plot +
    ggplot2::geom_point(data = dat_plot_rond, ggplot2::aes(x = outcome, y = exposure, size = log10_pval, color = z_score)) +
    ggplot2::scale_color_gradient2(name = name_zcore,
                                   low = scales::muted("#5884E5"),
                                   mid = "white",
                                   high = scales::muted("#9E131E"),
                                   midpoint = 0
    ) +
    ggplot2::scale_shape_manual(name = "", values = c(4,1)) +
    # ggplot2::scale_size(name = expression(-Log[10](P))) +
    
    ggplot2::scale_size(name = name_pval) +
    ggplot2::coord_fixed(clip = "off", ratio = 1) +
    ggplot2::guides(size = guide_legend(order = 1),
                    shape = guide_legend(order = 2)) +
    ggplot2::theme(
      panel.grid.major.y = element_line(size = 0.25, colour = "gray60"),
      panel.grid.major.x = element_line(size = 0.25, colour = "gray60"),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(t = 2, r = 0.5, b = 0.5, l = 0.5, "cm"),
      legend.position = "right",
      legend.text = element_text(
        color = "gray20",
        size = 10,
        margin = margin(l = 0.2, r = 0.2)
      ),
      legend.spacing.y = unit(0.1, 'cm'),
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      legend.key.size = unit(0.8, "cm"),
      axis.title = element_blank(),
      axis.line = element_line(linewidth = 0.5, colour = "gray20"),
      axis.ticks = element_line(linewidth = 0.5, colour = "gray20"),
      axis.text.y = element_text(
        size = 10,
        colour = "gray20"
      ),
      axis.text.x = element_text(
        angle = 60,
        size = 8,
        hjust = 1,
        face = "plain",
        colour = "gray20"
      ))
  
  if(subdivide_by_category){
    for (i in 1:length(ordered_category$Category[!is.na(ordered_category$Category)])) {
      
      balloon_plot = balloon_plot + ggplot2::annotation_custom(grob =  grid::textGrob(label = ordered_category$Category[!is.na(ordered_category$Category)][i],
                                                                                      hjust = 0,
                                                                                      vjust = 0.5,
                                                                                      rot = 18,
                                                                                      gp = grid::gpar(cex = 0.75, fontface = "bold")),
                                                               ymin = length(unique(dat_plot$exposure)) + 1,
                                                               ymax = length(unique(dat_plot$exposure)) + 1,
                                                               xmin = label_breaks[i],
                                                               xmax = label_breaks[i])
    }
  }
  
  print(balloon_plot)
  
}

  