library(data.table)
library(TAWApost)

###########################################################
# 
###########################################################

inflate_by_cpi <- function(tawa_data, inflator_file, base_year){
  
  # read data
  inflators <- fread(inflator_file)
  setnames(inflators, "year", "Year")
  
  # cpi for base year
  cpi_base <- inflators[Year == base_year, cpi]
  
  # cpi for each year
  tawa_data <- 
    merge(tawa_data, inflators[, .(Year, cpi)], by = "Year")
  
  # inflate
  tawa_data[, Adjust := cpi_base/cpi]
  tawa_data[, ':=' (HH_DI_BHC_Real = HH_DI_BHC*Adjust,
                    HH_DI_AHC_Real = HH_DI_AHC*Adjust,
                    HEDI_BHC_Real = HEDI_BHC*Adjust,
                    HEDI_AHC_Real = HEDI_AHC*Adjust)]
  
  # tidy 
  tawa_data[, ':=' (cpi=NULL, Adjust=NULL)]
  
  # better names 
  setnames(tawa_data, 
           c("HH_DI_BHC_Real",  "HH_DI_AHC_Real", "HEDI_BHC_Real", "HEDI_AHC_Real"), 
           paste0(c("HH_DI_BHC_",  "HH_DI_AHC_", "HEDI_BHC_", "HEDI_AHC_"), paste0(base_year, "dollars")))
  
  return(tawa_data)

}

###########################################################
# Distribution functions
###########################################################

# calculate precision based on stats minimum counts
calculate_precision <- function(n_samples){
  cut_points <- c(0, 10, 20, 50, 100, 500)
  min_n <- max(cut_points[n_samples>=cut_points])
  precision <- 5/max(5,min_n)
  return(precision)
}

# to calculate quantiles to highest precision 
calculate_quantiles <- function(tawa_data, value_name, group_name = "All"){
  
  # select group
  if (group_name ==  "All"){
    tawa_data[, Group := "All"]
  } else {
    tawa_data[, Group := .SD, .SDcols = group_name] 
  }
  
  # select value
  tawa_data[, Value := .SD, .SDcols = value_name] 
  
  # work out best precision
  tawa_data[, P_ID := .I] # not sure why i have to do this to get it to work? 
  tawa_data[, Samples := .N, by = .(Year, Group)]
  tawa_data[, Precision := calculate_precision(Samples), by = P_ID]
  
  # calculate weighted quantiles
  Quantiles <- 
    tawa_data[, .(Value = 
                    weighted_quantile(Value,Weight,
                                      seq(min(max(Precision),1-max(Precision)),
                                          max(max(Precision),1-max(Precision)), 
                                          max(Precision))), 
                  Quantile = seq(min(max(Precision),1-max(Precision)),
                                 max(max(Precision),1-max(Precision)), 
                                 max(Precision)), # need this strange workaround to make the function work when we dont have enough samples
                  Sample = max(Samples), 
                  Population = sum(Weight), 
                  Value_type = value_name,
                  Group_type = group_name), 
              by = .(Year, Group)]
  
  # remove the 0 and 1 percentiles (just there because of the work around with seq)
  Quantiles <- Quantiles[!(Quantile %in% c(0,1))]
  return(Quantiles)
}

convert_quantiles_to_histogram <- function(data_i, bins, min_Value, max_Value){
  
  if (data_i[!is.na(Value),.N]>2){
    
    data_i[,Cumulative_Population := Quantile*Population]
    
    total_population <- data_i[, max(Population)]
    
    if (length(bins) == 1){
    fit <- approx(c(min_Value,data_i[, Value], max_Value), 
                  c(0,data_i[, Cumulative_Population], total_population),
                  n=bins, 
                  ties="ordered")
    } else {
      fit <- approx(c(min_Value,data_i[, Value], max_Value), 
                    c(0,data_i[, Cumulative_Population], total_population),
                    xout=bins, 
                    ties="ordered")      
    }
    
    fit_data <- data.table(Vi=fit$x, 
                           Pi=fit$y)
    
    top_bin_population <- total_population-fit_data[, max(Pi)]
    
    fit_data[, ':='
             (Population=c(diff(Pi),top_bin_population))]
    
    return(fit_data[, .(Vi, Pi, Population)])
  }
}

convert_quantiles_to_mean <- function(quantile_data){
  
  quantile_data[, Midpoint := Value+c(diff(Value), Value[.N]-Value[.N-1])/2,by = .(Year, Group, Value_type, Group_type)]
  mean_data <- quantile_data[, .(Mean = sum(Midpoint*Quantile[1]), Quantile = Quantile[1]), by = .(Year, Group, Value_type, Group_type)]
  
  return(mean_data)
}

convert_quantiles_to_population_below_threshold <- function(quantile_data, thresholds){
  
  quantile_data[, 
                Cumulative_Population := Quantile*Population]
  
  below_threshold_data <- data.table()
  
  for (y in quantile_data[, unique(Year)]){
    for (g in quantile_data[, unique(Group)]){
      data_i <- quantile_data[Year == y &  Group == g ]
      
      if (data_i[!is.na(Value),.N]>2){     
        total_population <- quantile_data[Year == y & Group == g, max(Population)]
        max_Value <- 10^10
        threshold <- thresholds[Year==y, Threshold] 
        
        fit <- approx(c(0,quantile_data[Year == y & Group == g,Value], max_Value), 
                      c(0,quantile_data[Year == y & Group == g, Cumulative_Population], total_population),
                      threshold, 
                      ties="ordered")
        
        fit_data <- data.table(Threshold=fit$x, 
                               Population=fit$y)
        
        fit_data[, ':='
                 (Proportion=Population/total_population, 
                   Year=y, 
                   Group=g)]
        
        below_threshold_data <- rbind(below_threshold_data, 
                                      fit_data[, .(Threshold, Population, Proportion, Year, Group)])
      }
    }
  }
  
  return(below_threshold_data)
}

sevenyear_comparisons <- function(input_data, type){
  
  if (type == "Year"){
    Change <- input_data[ Year %in% c(2007, 2014, 2021)]
    
    Change <- merge(Change,
                    Change[Year == 2007, .( 
                      Population_2007 = Population, 
                      Group_Total_2007 = Group_Total, 
                      All_Total_2007 = All_Total, 
                      Income, 
                      Group)], 
                    by = c("Income", "Group"), 
                    allow.cartesian = TRUE) 
    
    Change <- merge(Change,
                    Change[Year == 2014, .( 
                      Population_2014 = Population, 
                      Group_Total_2014 = Group_Total, 
                      All_Total_2014 = All_Total, 
                      Income, 
                      Group)], 
                    by = c("Income", "Group"), 
                    allow.cartesian = TRUE) 
    
    Change[, Diff_prop_population := Population/All_Total - Population_2007/All_Total_2007]
    Change[, Diff_prop_group := Population/Group_Total - Population_2007/Group_Total_2007]
    
    Change[Year == 2021, Diff_prop_population := Population/All_Total - Population_2014/All_Total_2014]
    Change[Year == 2021, Diff_prop_group := Population/Group_Total - Population_2014/Group_Total_2014]
  } 
  
  if(type == "Year_group"){
    Change <- input_data[Year_group %in% c("2007-2009", "2013-2015", "2021")]
    
    Change <- merge(Change,
                    Change[Year_group == "2007-2009", .( 
                      Population_2007 = Population, 
                      Group_Total_2007 = Group_Total, 
                      All_Total_2007 = All_Total, 
                      Income, 
                      Group)], 
                    by = c("Income", "Group"), 
                    allow.cartesian = TRUE) 
    
    Change <- merge(Change,
                    Change[Year_group == "2013-2015", .( 
                      Population_2014 = Population, 
                      Group_Total_2014 = Group_Total, 
                      All_Total_2014 = All_Total, 
                      Income, 
                      Group)], 
                    by = c("Income", "Group"), 
                    allow.cartesian = TRUE) 
    
    Change[, Diff_prop_population := Population/All_Total - Population_2007/All_Total_2007]
    Change[, Diff_prop_group := Population/Group_Total - Population_2007/Group_Total_2007]
    
    Change[Year_group == "2021", Diff_prop_population := Population/All_Total - Population_2014/All_Total_2014]
    Change[Year_group == "2021", Diff_prop_group := Population/Group_Total - Population_2014/Group_Total_2014]
  }
    
  return(Change)
  
}


plot_inequality_metrics <- function(decile_data, gini_data = NA){
  
  plot_inequality <- data.table()
  
  if(!is.na(gini_data)){
    plot_inequality <- rbind(plot_inequality, 
                             gini_data)
  }
  
  ratios <- dcast(decile_data, Year + Income_type + Group ~ Statistic, value.var = "Value")
  setnames(ratios , gsub(" ", "_", names(ratios)))
  
  ratios[, ':=' (R_9010 = Decile_9/Decile_1, 
                 R_8020 = Decile_8/Decile_2)]
  ratios <- ratios[,.(Year, Group, R_9010, R_8020, Income_type)]
  ratios <- melt(ratios,
                 id.vars = c("Year", "Income_type", "Group"),
                 variable.name = "Statistic", 
                 value.name = "Value")

  plot_inequality <- rbind(plot_inequality, 
                           ratios, 
                           fill=T)
  
  plot_inequality[Income_type == "BHC", Income_type := "Before housing cost"]
  plot_inequality[Income_type == "AHC", Income_type := "After housing cost"]
  plot_inequality[Statistic == "R_9010", Statistic := "90/10"]
  plot_inequality[Statistic == "R_8020", Statistic := "80/20"]
  
  p <- 
    ggplot(plot_inequality,
         aes(x=Year,
             y=Value, 
             color=Group)) +
    geom_line() + 
    facet_grid(Statistic~Income_type,scales="free") + 
    theme_bw() +   
    theme(strip.background = element_rect(colour="#000000", fill="#FFFFFF")) + 
    scale_y_continuous(name="") 
    

  p <- ggplotly(p)
  return(p)
}


plot_medians <- function(decile_data){
  plot_levels <- copy(decile_data)
  plot_levels <- plot_levels[Statistic == "Decile 5"]
  
  plot_levels[Income_type == "BHC", Income_type := "Before housing cost"]
  plot_levels[Income_type == "AHC", Income_type := "After housing cost"]
  
  p <- ggplot(plot_levels,
         aes(x=Year,
             y=Value,
             color=Group)) +
    geom_line() + 
    facet_grid(.~Income_type,scales="free") + 
    theme_bw() +   
    theme(strip.background = element_rect(colour="#000000", fill="#FFFFFF")) + 
    scale_y_continuous(labels=dollar_format(), 
                       name="Median")
  
  p <- ggplotly(p)
  
  return(p)
  
}


average_annual_change_deciles <- function(decile_data){
  decile_lag <- copy(decile_data)
  decile_lag[, Year := Year+1]
  decile_lag[, Previous_Value := Value]
  decile_lag[, Value := NULL]
  
  decile_all <- merge(decile_data, 
                      decile_lag, 
                      by = c("Year", "Statistic", "Income_type", "Group"))
  decile_all[, Prop_Increase := Value/Previous_Value-1]
  decile_all[, Dollar_Increase := Value-Previous_Value]
  
  decile_all <- 
    decile_all[, .(Prop_Increase = exp(mean(log(1+Prop_Increase)))-1, 
                   Dollar_Increase = mean(Dollar_Increase)), 
               by = .(Statistic, Income_type, Group)]
  
  decile_all <- melt(decile_all, id.vars = c("Statistic", "Income_type", "Group"))
  decile_all <- dcast(decile_all, 
                      Group + Statistic ~ Income_type + variable)
  
  decile_all[, ':=' (AHC_Prop_Increase=percent(AHC_Prop_Increase, accuracy=.01), 
                     BHC_Prop_Increase=percent(BHC_Prop_Increase, accuracy=.01), 
                     AHC_Dollar_Increase=dollar(AHC_Dollar_Increase, accuracy=1), 
                     BHC_Dollar_Increase=dollar(BHC_Dollar_Increase, accuracy=1))]
  
  
  decile_table <- 
    decile_all %>% 
    kbl(col.names = c(" ","", "Proportion", "Dollar", "Proportion", "Dollar")) %>% 
    kable_styling() %>% 
    add_header_above(c(" " = 2, 
                       "Average annual AHC increase" = 2, 
                       "Average annual BHC increase" = 2))  %>%
    collapse_rows(columns = 1, valign = "top")
  
  return(decile_table)
}



###########################################################
# Groupings
###########################################################

#' assign a factor to each individual in a group 
#' with values <18, 18-30, 31-40, 41-50, 51-64, 65+
#'
#' @param age vector of ages of people in group
#' 
#' @return 
#' @export
#'
#' @examples
age_group <- function(age){
  this_age_group <- cut(age, 
                        breaks = c(-1,17,30,40,50,64,200), 
                        labels = FALSE)
  this_age_group <- factor(this_age_group, 
                           c(1,2,3,4,5,6), 
                           labels = c("<18", "18-30", "31-40", "41-50", "51-64", "65+"))
  return(this_age_group)
}

#' age group of the youngest child (where a child is defined as under 18 years old) in a group 
#' with values <3, 3-4, 5-11, 12+, No children
#'
#' @param age vector of ages of people in group
#' 
#' @return 
#' @export
#'
#' @examples
youngest_child_in_group <- function(age){
  youngest_child <- cut(min(age), 
                        breaks = c(-1,2,4,11,17,200), 
                        labels = FALSE)
  youngest_child <- factor(youngest_child, 
                           c(1,2,3,4,5), 
                           labels=c("<3", "3-4", "5-11", "12+", "No children"))
  return(youngest_child)
}

#' age group of the oldest person in a group 
#' with values <=30, 31-40, 41-50, 51-64, 65+
#'
#' @param age vector of ages of people in group
#' 
#' @return 
#' @export
#'
#' @examples
oldest_in_group <- function(age){
  oldest_person <- cut(max(age), 
                       breaks = c(-1,30,40,50,64,200), 
                       labels = FALSE)
  oldest_person <- factor(oldest_person, 
                          c(1,2,3,4,5), 
                          labels = c("<=30", "31-40", "41-50", "51-64", "65+"))
  return(oldest_person)
}

#' birth decade of oldest 
#' with values <1950s, 1950-1959, 1960-1969,  1970-1979,  1980-1989,  1990-1999,  2000-2009,  2010-2019
#'
#' @param age vector of birth years of people in group
#' 
#' @return 
#' @export
#'
#' @examples
decade_of_oldest <- function(birth_year){
  oldest_decade <- cut(min(birth_year), 
                       breaks = c(-1,seq(1950, 2020, 10),10000), 
                       labels = FALSE,
                       right=FALSE)
  oldest_decade <- factor(oldest_decade, 
                          1:9, 
                          labels = c("<1950s", "1950-1959", "1960-1969",  "1970-1979",  "1980-1989",  "1990-1999",  "2000-2009",  "2010-2019", "2020+"))
  return(oldest_decade)
}


#' calculates household type of a group based on an age vector and 
#' family id
#' 
#' @param age vector of ages
#' @param family_id
#' 
#' @return 
#' @export
#'
#' @examples
#' 
household_type <- function(age, family_id){
  
  fams <- length(unique(family_id))
  adults <- sum(age>=18)
  children <- sum(age<18)
  super <- sum(age>=65)
  
  if (fams>1){
    if(children==0){
      hh_type <- "Multi-family"
    } else {
      hh_type <- "Multi-family with children"
    }
  } else {
    if(adults==1){
      if(children==0){
        hh_type <- "Single"
      } else {
        hh_type <- "Single parent"
      }
    } else {
      if(children==0){
        hh_type <- "Couple"
      } else {
        hh_type <- "Couple parents"
      }  
    }
  }
  
  if (super>0){
    if(adults==1){
      hh_type <- "65+, single"
      } else {
      hh_type <- "65+, couple/sharing"
      }
  } 
  hh_type <- factor(hh_type,
                    levels = c( "Single","Single parent",
                                "Couple","Couple parents",
                                "Multi-family","Multi-family with children",
                                "65+, single","65+, couple/sharing"))
  return(hh_type)
}

