## Dernière date de mise à jour : 02/08/2024 ##

#### ggcompetingrisks1 ####
ggcompetingrisks1 <- function(fit, gnames = NULL, gsep=" ", 
                              multiple_panels = FALSE,
                              event_suppr=NULL, 
                              lwd=0.5, 
                              palette=NULL, 
                              labs=NULL, labs_event=NULL, 
                              type_group="color", type_event="linetype", linetype=NULL,
                              alpha = 0.05, conf.int = FALSE, 
                              ggtheme = theme_survminer(), ...) {
  
  
  
  #' Cumulative Incidence Curves for Competing Risks
  #' @description This function plots Cumulative Incidence Curves. For \code{cuminc} objects it's a \code{ggplot2} version of \code{plot.cuminc}.
  #' @param fit an object of a class \code{cmprsk::cuminc} - created with \code{cmprsk::cuminc} function.
  #'   Warning, the event variable must be without spaces.
  #' @param gnames a vector with group names. If not supplied then will be extracted from \code{fit} object (\code{cuminc} only).
  #' @param gsep a separator that extracts group names and event names from \code{gnames} object (\code{cuminc} only).
  #' @param multiple_panels if \code{TRUE} then groups will be plotted in different panels (\code{cuminc} only).
  #' @param event_suppr value of events that you don't want to plot
  #' @param lwd line thickness
  #' @param palette the color palette to be used
  #'  Default is palette from hue_pal
  #' @param labs character vector specifying legend labels. Used to replace the names of the strata from the fit. Should be given in the same order as those strata.
  #' @param labs_event character vector specifying legend labels for the events if more than one. Used to replace the names of the events from the fit. Should be given in the same order as those events
  #' @param type_group if there is more than one group, used to specify if you want to differentiate the groups by different colors ("color"), different line types ("linetype") or both ("color_linetype"). 
  #'  Default is \code{"color"}
  #' @param type_event if there is more than one event, used to specify if you want to differentiate the events by different colors ("color"), different line types ("linetype") or both ("color_linetype"). 
  #'  Default is \code{"linetype"}
  #'  If more than one group and more than one event, "color" will be used for \code{type_group} and "linetype" will be used for \code{type_event}.
  #' @param alpha see \code{conf.int}, scaling actor for the ribbon. The default value is 0.05.
  #' @param conf.int if \code{TRUE} then additional layer (\code{geom_ribbon}) is added around the point estimate. The ribon is plotted with boundries +- \code{qnorm(1-alpha/2)}*standard deviation. 
  #'  standard deviation is estimated with cif.ci method
  #'  Default is \code{FALSE}
  #' @param ggtheme function, \code{ggplot2} theme name. Default value is \link{theme_survminer}.
  #'  Allowed values include ggplot2 official themes: see \code{\link[ggplot2]{theme}}.
  #' @param legend character specifying legend position. Allowed values are one of c("top", "bottom", "left", "right", "none"). 
  #'  Default is "top" side position. to remove the legend use legend = "none". Legend position can be also specified using a numeric vector c(x, y). 
  #'  In this case it is possible to position the legend inside the plotting area. x and y are the coordinates of the legend box. 
  #'  Their values should be between 0 and 1. c(0,0) corresponds to the "bottom left" and c(1,1) corresponds to the "top right" position. For instance use legend = c(0.8, 0.2).
  #' @param legend.title legend title.
  #' @param ... further arguments passed to the function \code{\link[ggpubr]{ggpar}} for customizing the plot.
  #' @return Returns an object of class \code{gg}.
  #'
  #'
  #'
  
  if(any(!class(fit) %in% c("cuminc"))){
    stop("Erreur : l'object n'est pas de la classe cuminc")
  }

  if (any(class(fit) == "cuminc")) {
    pl <- ggcompetingrisks.cuminc1(fit = fit, gnames=gnames, gsep=gsep, 
                                   multiple_panels=multiple_panels, event_suppr=event_suppr, 
                                   lwd=lwd, palette=palette, 
                                   labs=labs, labs_event=labs_event,
                                   type_group=type_group, type_event=type_event, linetype=linetype,
                                   alpha = alpha, conf.int = conf.int)
  }
  
  pl <- pl + ggtheme +
    ylab("Cumulative incidence") + xlab("Time") +
    ggtitle("")
  ggpubr::ggpar(pl, ...)
}


ggcompetingrisks.cuminc1 <- function(fit, gnames = NULL, gsep=" ", 
                                     multiple_panels = FALSE, 
                                     event_suppr=NULL, 
                                     lwd=0.5,
                                     palette=NULL, 
                                     labs=NULL, labs_event=NULL,
                                     type_group="color", type_event="linetype", linetype=NULL,
                                     alpha = 0.05, conf.int = FALSE) {
  if (!is.null(fit$Tests))
    fit <- fit[names(fit) != "Tests"]
  fit2 <- lapply(fit, `[`, 1:3)
  if (is.null(gnames)) gnames <- names(fit2)
  fit2_list <- lapply(seq_along(gnames), function(ind) {
    df <- as.data.frame(fit2[[ind]])
    df$name <- gnames[ind]
    df
  })

  df <- do.call(rbind, fit2_list)
  
  require(msm)
  
  df$lse <- apply(cbind(df$est,df$var),1,function(x){deltamethod(~log(-log(x1)),x[1],x[2])})
  
  df$lse[is.na(df$lse)] <- 0

  df$low <- exp(-exp(log(-log(df$est))+qnorm(1-alpha/2)*df$lse))
  df$upp <- exp(-exp(log(-log(df$est))-qnorm(1-alpha/2)*df$lse))
  
  
  
  
  
  # # probleme si var groupe avec des espaces
  # df$event <- sapply(strsplit(df$name, split=gsep), `[`, 2)
  # df$group <- sapply(strsplit(df$name, split=gsep), `[`, 1)
  
  
  tmp <- data.frame(stringr::str_split(df$name, gsep, simplify = T))
  df$event <- tmp[, ncol(tmp)]
  
  if(ncol(tmp)==2){
    df$group <- tmp[, -ncol(tmp)]
  }else{
    df$group <- apply(tmp[, -ncol(tmp)], 1, paste, collapse=" ")
    warning(paste0("These events have been taken into account: ", paste(unique(df$event), collapse = ", "), ". If incorrect, please recode the event without spaces."))
  }

  
  
  df$std <- std <- sqrt(df$var)
  
  df$group <- factor(df$group, unique(df$group))
  
  # we delete the events we don't want to plot
  if(!is.null(event_suppr)){
    df=df[!df$event %in% event_suppr, ]
    if(length(unique(df$event))==1) {df$event=NULL}
  }
  
  
  if(length(unique(df$group))>1 & !is.null(df$event)){
    warning(paste0("More than one group and more than one event to plot: groups have been differentiated by color and events by line type"))
    type_group <- "color"
    type_event <- "linetype"
  }

  
  
  # if palette argument is empty
  if(is.null(palette) & type_group %in% c("color", "color_linetype") & length(unique(df$group))!=1){
    palette <- scales::hue_pal()(length(unique(df$group)))
  }
  if(is.null(palette) & type_event %in% c("color", "color_linetype") & !is.null(df$event)){
    palette <- scales::hue_pal()(length(unique(df$event)))
  }
  
  # if palette argument does not contain enough elements
  if(!is.null(palette) & type_group %in% c("color", "color_linetype") & length(unique(df$group))!=1 & length(palette)<length(unique(df$group)) & length(palette)){
    if(!(length(palette)==1 & palette[1]=="hue")){
      warning(paste0("The palette argument does not contain enough elements: ", length(unique(df$group)) , " needed but only ", length(palette) , " provided."))
    }
    palette <- scales::hue_pal()(length(unique(df$group)))
  }
  if(!is.null(palette) & type_event %in% c("color", "color_linetype") & length(unique(df$group))==1 & length(unique(df$event))!=1 & length(palette)<length(unique(df$event))){
    if(!(length(palette)==1 & palette[1]=="hue")){
      warning(paste0("The palette argument does not contain enough elements: ", length(unique(df$event)) , " needed but only ", length(palette) , " provided. hue_pal() has been used"))
    }
    palette <- scales::hue_pal()(length(unique(df$event)))
  }
  
  # if palette argument contain too many elements
  if(!is.null(palette) & type_group %in% c("color", "color_linetype") & length(unique(df$group))!=1 & length(palette)>length(unique(df$group))){
    palette <- palette[1:length(unique(df$group))]
    warning(paste0("The palette argument contains too many elements: only the first elements of the vector have been used (", palette, ")"))
  }
  if(!is.null(palette) & type_event %in% c("color", "color_linetype") & length(unique(df$group))==1 & length(unique(df$event))!=1 & length(palette)>length(unique(df$event))){
    palette <- palette[1:length(unique(df$event))]
    warning(paste0("The palette argument contains too many elements: only the first elements of the vector have been used (", palette, ")"))
  } 
  
  
  # if linetype argument is empty
  if(is.null(linetype) & type_group %in% c("linetype", "color_linetype") & length(unique(df$group))!=1){
    linetype <- 1:length(unique(df$group))
  }
  if(is.null(linetype) & type_event %in% c("linetype", "color_linetype") & !is.null(df$event)){
    linetype <- 1:length(unique(df$event))
  }
  
  # if linetype argument does not contain enough elements
  if(!is.null(linetype) & type_group %in% c("linetype", "color_linetype") & length(unique(df$group))!=1 & length(linetype)<length(unique(df$group))){
    warning(paste0("The linetype argument does not contain enough elements: ", length(unique(df$group)) , " needed but only ", length(linetype) , " provided."))
    linetype <- 1:length(unique(df$group))
  }
  if(!is.null(linetype) & type_event %in% c("linetype", "color_linetype") & length(unique(df$group))==1 & length(unique(df$event))!=1 & length(linetype)<length(unique(df$event))){
    warning(paste0("The linetype argument does not contain enough elements: ", length(unique(df$event)) , " needed but only ", length(linetype) , " provided. hue_pal() has been used"))
    linetype <- 1:length(unique(df$event))
  }
  
  # if linetype argument contain too many elements
  if(!is.null(linetype) & type_group %in% c("linetype", "color_linetype") & length(unique(df$group))!=1 & length(linetype)>length(unique(df$group))){
    linetype <- linetype[1:length(unique(df$group))]
    warning(paste0("The linetype argument contains too many elements: only the first elements of the vector have been used (", linetype, ")"))
  }
  if(!is.null(linetype) & type_event %in% c("linetype", "color_linetype") & length(unique(df$group))==1 & length(unique(df$event))!=1 & length(linetype)>length(unique(df$event))){
    linetype <- linetype[1:length(unique(df$event))]
    warning(paste0("The linetype argument contains too many elements: only the first elements of the vector have been used (", linetype, ")"))
  } 
  
  
  
  
  
  

  
  
  
  if(length(unique(df$group))==1){
    
    if(is.null(df$event)) {
      
      if(length(palette)>1){
        palette <- palette[1]
        warning(paste0("The palette argument contains more than one element : only the first element of the vector has been used (", palette, ")"))
      }
        
      
      pl <- ggplot(df, aes(time, est, color=palette), legend="none")
      
      if(conf.int) pl <- pl + geom_ribbon(aes(ymin = low, ymax=upp, fill=palette), alpha = 0.3, linetype=0)
      
    }else{
      
      if(type_event=="color"){
        pl <- ggplot(df, aes(time, est, color=event))
        if(conf.int) pl <- pl + geom_ribbon(aes(ymin = low, ymax=upp, fill=event), alpha = 0.3, linetype=0) 
      }
      if(type_event=="linetype"){
        if(is.null(palette)){
          palette <- "black"
        }
        if(length(palette)==1 & palette[1]=="hue"){
          palette <- scales::hue_pal()(1)
        }
        pl <- ggplot(df, aes(time, est, linetype=event, color=event))
        if(conf.int) {
          pl <- pl + geom_ribbon(aes(ymin = low, ymax=upp, fill=event), alpha = 0.3, linetype=0) 
        }
      }
      if(type_event=="color_linetype"){
        pl <- ggplot(df, aes(time, est, linetype=event, color=event))
        if(conf.int) {
          pl <- pl + geom_ribbon(aes(ymin = low, ymax=upp, fill=event), alpha = 0.3, linetype=0) 
        }
      }

      }
    
    
  } else {
    if (multiple_panels) {
      pl <- ggplot(df, aes(time, est, color=event)) + facet_wrap(~group)
      if(conf.int) pl <- pl + geom_ribbon(aes(ymin = low, ymax=upp, fill=palette), alpha = 0.3, linetype=0) 
    } else {
      if(is.null(df$event)){
        if(type_group=="color"){
          pl <- ggplot(df, aes(time, est, color=group))
          if(conf.int) pl <- pl + geom_ribbon(aes(ymin = low, ymax=upp, fill=group), alpha = 0.3, linetype=0) 
        }
        if(type_group=="linetype"){
          if(is.null(palette)){
            palette <- "black"
          }
          if(length(palette)==1 & palette[1]=="hue"){
            palette <- scales::hue_pal()(1)
          }
          pl <- ggplot(df, aes(time, est, linetype=group, color=group))
          if(conf.int) {
            pl <- pl + geom_ribbon(aes(ymin = low, ymax=upp, fill=group), alpha = 0.3, linetype=0) 
          }
        }
        if(type_group=="color_linetype"){
          pl <- ggplot(df, aes(time, est, linetype=group, color=group))
          if(conf.int) {
            pl <- pl + geom_ribbon(aes(ymin = low, ymax=upp, fill=group), alpha = 0.3, linetype=0) 
          }
        }
      } else {
        pl <- ggplot(df, aes(time, est, color=group, linetype=event))
        if(conf.int){
          for(i in unique(df$event)){
            pl <- pl + geom_ribbon(data = df[df$event==i, ], aes(ymin = low, ymax=upp, fill=group), alpha = 0.3, linetype=0) 
          }

        }
      }

      
    }
  }
  
  
  
  
  
  if(is.null(labs)){labs <- unique(df$group)}
  if(is.null(labs_event)){labs_event <- unique(df$event)}
  
  
  if(is.null(df$event)){
    if(type_group!="linetype"){
      pl <- pl +
        scale_color_manual(values = palette, labels = labs) +
        scale_fill_manual(values = palette, labels = labs) +
        scale_linetype_manual(values = linetype, labels = labs)
    } else{
      pl <- pl +
        scale_color_manual(values = rep(palette, length(unique(df$group))), labels = labs)+
        scale_fill_manual(values = rep(palette, length(unique(df$group))), labels = labs) +
        scale_linetype_manual(values = linetype, labels = labs)
    }
  } else if (length(unique(df$group))>1){
    pl <- pl +
      scale_color_manual(values = palette, labels = labs) +
      scale_fill_manual(values = palette, labels = labs) +
      scale_linetype_manual(values = linetype, labels = labs_event)
  } else{
    if(type_event!="linetype"){
      pl <- pl +
        scale_color_manual(values = palette, labels = labs_event) +
        scale_fill_manual(values = palette, labels = labs_event) +
        scale_linetype_manual(values = linetype, labels = labs_event)
    } else{
      pl <- pl +
        scale_color_manual(values = rep(palette, length(unique(df$event))), labels = labs_event)+
        scale_fill_manual(values = rep(palette, length(unique(df$event))), labels = labs_event) +
        scale_linetype_manual(values = linetype, labels = labs_event)
    }
  }


  
  
  pl + geom_line(lwd=lwd) 
 
}







#### format.pv ####
format.pv <- function(p, text=F)
{
  if(p<0.0001) return("<0.0001")
  if(p>=0.0001&p<0.00095) ifelse(text==F,return(sprintf("%.4f", p)),return(paste("=",sprintf("%.4f", p),sep="")))
  if(p>=0.00095&p<=0.0095) ifelse(text==F,return(as.character(signif(p,1))),return(paste("=",as.character(signif(p,1)),sep="")))
  if(p>0.0095&p<0.0995) ifelse(text==F,return(sprintf("%.3f", signif(p,2))),return(paste("=",sprintf("%.3f", signif(p,2)),sep="")))
  if(p>=0.0995) ifelse(text==F,return(sprintf("%.2f", signif(p,2))),return(paste("=",sprintf("%.2f", signif(p,2)),sep="")))
}





#### ggcombine ####
ggcombine <- function(list_obj, name_obj, linetype_obj=NULL){
  
  #' Combine ggsurvplot and ggcompeting risks
  #' @description This function plots ggsurvplot and ggcompetingrisks in a same graphic
  #' @param list_obj list of objects of a class ggsurvplot or ggcompetingrisks1
  #' @param name_obj a vector with plot names
  #' @param linetype_obj a vector with plot groups (to differenciate groups by linetype, only for cumulative incidence)
  #' @return Returns an object of class \code{gg}.
  #'
  #'
  
  
  plot.all <- ggplot() 
  
  vect_color <- NULL
  
  for(i in 1:length(list_obj)){
    
    if(length(grep("surv", names(list_obj[[i]])))>0){
      
      data_surv <- ggplot_build(list_obj[[i]]$plot)$data[[1]]
      data_surv$name <- name_obj[i]
      
      if(length(ggplot_build(list_obj[[i]]$plot)$data)==3){ # without confint
        
        data_cens <- ggplot_build(list_obj[[i]]$plot)$data[[3]]
        data_cens$name <- name_obj[i]
        
      }else if(length(ggplot_build(list_obj[[i]]$plot)$data)==4){ # with confint
        
        data_confint <- ggplot_build(list_obj[[i]]$plot)$data[[3]]
        data_confint$name <- name_obj[i]
        data_cens <- ggplot_build(list_obj[[i]]$plot)$data[[4]]
        data_cens$name <- name_obj[i]
        
      }
      
      
      id_color <- unique(data_surv$colour)
      
      plot.all <- plot.all +
        geom_step(data = data_surv, aes(x, y, color=name)) + 
        geom_point(data = data_cens, aes(x, y, color=name), shape=3)
      
      
    }else{
      data_icc <- ggplot_build(list_obj[[i]])$data[[1]]
      data_icc$name <- name_obj[i]
      
      id_color <- unique(data_icc$colour)
      
      
      if(!is.null(linetype_obj)){
        data_icc$linetype <- linetype_obj[i]
        plot.all <- plot.all +
          geom_step(data = data_icc, aes(x, y, color=name, linetype=linetype))
      }else{
        plot.all <- plot.all +
          geom_step(data = data_icc, aes(x, y, color=name))
      }
      
    }
    
    # print(plot.all)
    
    
    vect_color <- c(vect_color, setNames(id_color, name_obj[i]))
    
    
  }
  
  plot.all <- plot.all +
    scale_color_manual(values=vect_color)
  
  return(plot.all)
  
}

