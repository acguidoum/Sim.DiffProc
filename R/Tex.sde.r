## Thu Apr 30 10:50:58 2020
## Original file Copyright © 2020 A.C. Guidoum, K. Boukhetala
## This file is part of the R package Sim.DiffProc
## Department of Probabilities & Statistics
## Faculty of Mathematics
## University of Science and Technology Houari Boumediene
## BP 32 El-Alia, U.S.T.H.B, Algiers
## Algeria

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## A copy of the GNU General Public License is available at
## http://www.r-project.org/Licenses/
## Unlimited use and distribution (see LICENCE).
###################################################################################################


#####
##### LaTeX code for Sim.DiffProc package

TEX.sde <- function(object, ...)  UseMethod("TEX.sde")

TEX.sde.default <- function(object, ...)
        {

   greek <- c("alpha", "theta", "tau", "beta", "vartheta", "pi", "upsilon",
           "gamma", "varpi", "phi", "delta", "kappa", "rho","iota",
           "varphi", "epsilon", "lambda", "varrho", "chi", "varepsilon",
           "mu", "sigma", "psi", "zeta", "nu", "varsigma", "omega", "eta",
           "xi", "Gamma", "Lambda", "Sigma", "Psi", "Delta", "Xi", 
           "Upsilon", "Omega", "Theta", "Pi", "Phi")
   greek0 <- c(paste0(greek,"0"))
   greek1 <- c(paste0(greek,"1"))
   greek2 <- c(paste0(greek,"2"))
   greek3 <- c(paste0(greek,"3"))
   greek4 <- c(paste0(greek,"4"))
   greek5 <- c(paste0(greek,"5"))
   greek6 <- c(paste0(greek,"6"))
   greek7 <- c(paste0(greek,"7"))
   greek8 <- c(paste0(greek,"8"))
   greek9 <- c(paste0(greek,"9"))
   greek10 <- c(paste0(greek,"10"))
   greek11 <- c(paste0(greek,"11"))
   greek12 <- c(paste0(greek,"12"))
   greek13 <- c(paste0(greek,"13"))
   greek14 <- c(paste0(greek,"14"))
   greek15 <- c(paste0(greek,"15"))
   greek16 <- c(paste0(greek,"16"))
   greek17 <- c(paste0(greek,"17"))
   greek18 <- c(paste0(greek,"18"))
   greek19 <- c(paste0(greek,"19"))
   greek20 <- c(paste0(greek,"20"))
   greek21 <- c(paste0(greek,"21"))
   greek22 <- c(paste0(greek,"22"))
   greek23 <- c(paste0(greek,"23"))
   greek24 <- c(paste0(greek,"24"))
   greek25 <- c(paste0(greek,"25"))
   greek26 <- c(paste0(greek,"26"))
   greek27 <- c(paste0(greek,"27"))
   greek28 <- c(paste0(greek,"28"))
   greek29 <- c(paste0(greek,"29"))
   greek30 <- c(paste0(greek,"30"))
   greek_list  <- setNames(paste0("\\", greek),greek)
   greek_list0 <- setNames(paste0("\\", greek,"_","{0}"), greek0)
   greek_list1 <- setNames(paste0("\\", greek,"_","{1}"), greek1)
   greek_list2 <- setNames(paste0("\\", greek,"_","{2}"), greek2)
   greek_list3 <- setNames(paste0("\\", greek,"_","{3}"), greek3)
   greek_list4 <- setNames(paste0("\\", greek,"_","{4}"), greek4)
   greek_list5 <- setNames(paste0("\\", greek,"_","{5}"), greek5)
   greek_list6 <- setNames(paste0("\\", greek,"_","{6}"), greek6)
   greek_list7 <- setNames(paste0("\\", greek,"_","{7}"), greek7)
   greek_list8 <- setNames(paste0("\\", greek,"_","{8}"), greek8)
   greek_list9 <- setNames(paste0("\\", greek,"_","{9}"), greek9)
   greek_list10 <- setNames(paste0("\\", greek,"_","{10}"), greek10)
   greek_list11 <- setNames(paste0("\\", greek,"_","{11}"), greek11)
   greek_list12 <- setNames(paste0("\\", greek,"_","{12}"), greek12)
   greek_list13 <- setNames(paste0("\\", greek,"_","{13}"), greek13)
   greek_list14 <- setNames(paste0("\\", greek,"_","{14}"), greek14)
   greek_list15 <- setNames(paste0("\\", greek,"_","{15}"), greek15)
   greek_list16 <- setNames(paste0("\\", greek,"_","{16}"), greek16)
   greek_list17 <- setNames(paste0("\\", greek,"_","{17}"), greek17)
   greek_list18 <- setNames(paste0("\\", greek,"_","{18}"), greek18)
   greek_list19 <- setNames(paste0("\\", greek,"_","{19}"), greek19)
   greek_list20 <- setNames(paste0("\\", greek,"_","{20}"), greek20)
   greek_list21 <- setNames(paste0("\\", greek,"_","{21}"), greek21)
   greek_list22 <- setNames(paste0("\\", greek,"_","{22}"), greek22)
   greek_list23 <- setNames(paste0("\\", greek,"_","{23}"), greek23)
   greek_list24 <- setNames(paste0("\\", greek,"_","{24}"), greek24)
   greek_list25 <- setNames(paste0("\\", greek,"_","{25}"), greek25)
   greek_list26 <- setNames(paste0("\\", greek,"_","{26}"), greek26)
   greek_list27 <- setNames(paste0("\\", greek,"_","{27}"), greek27)
   greek_list28 <- setNames(paste0("\\", greek,"_","{28}"), greek28)
   greek_list29 <- setNames(paste0("\\", greek,"_","{29}"), greek29)
   greek_list30 <- setNames(paste0("\\", greek,"_","{30}"), greek30)
   mem_sde <- c("m","S","m1","m2","m3","S1","S2","S3","C12","C13","C23")
   mem_sde_tex <- c("m(t)","S(t)","m_{1}(t)","m_{2}(t)","m_{3}(t)","S_{1}(t)","S_{2}(t)","S_{3}(t)","C_{12}(t)","C_{13}(t)","C_{23}(t)")
   mem_sde_list <- setNames(mem_sde_tex, mem_sde)
   var_sde <- c("x","y","z","w","w1","w2","w3","X","Y","Z","W","W1","W2","W3")
   var_sde_tex <- c("X_{t}","Y_{t}","Z_{t}","W_{t}","W_{1,t}","W_{2,t}",
                 "W_{3,t}","X_{t}","Y_{t}","Z_{t}","W_{t}","W_{1,t}","W_{2,t}",
                 "W_{3,t}")
   var_sde_list <- setNames(var_sde_tex, var_sde)
   #fct_sde <- c("P")
   #fct_sde_tex <- c("\\mathbb{P}")
   #fct_sde_list <- setNames(fct_sde_tex, fct_sde)
   greek_env   <- list2env(as.list(c(greek_list,greek_list0,greek_list1,greek_list2,greek_list3,greek_list4,greek_list5,
                                     greek_list6,greek_list7,greek_list8,greek_list9,greek_list10,greek_list11,greek_list12,
									 greek_list13,greek_list14,greek_list15,greek_list16,greek_list17,greek_list18,greek_list19,
									 greek_list20,greek_list21,greek_list22,greek_list23,greek_list24,greek_list25,greek_list26,
									 greek_list27,greek_list28,greek_list29,greek_list30,
                                     mem_sde_list,var_sde_list)), parent = emptyenv())
   unary_op <- function(left, right) {
           force(left)
           force(right)
           function(e1) {
              paste0(left, e1, right)
             }
     }
   # binary_op <- function(sep) {
           # force(sep)
           # function(e1, e2) {
              # paste0(e1, sep, e2)
             # }
     # }
   binary_op <- function(sep) {
           force(sep)
           function(e1, e2) {
            if(missing(e2)){ paste0(" - ", e1) }else{ paste0(e1,sep,e2)}
            }
     }
	 
   # Binary operators
   f_env <- new.env(parent = emptyenv())
   f_env$"+" <- binary_op(" + ")
   f_env$"-" <- binary_op(" - ")
   f_env$"*" <- binary_op(" \\, ")
   f_env$"/" <- binary_op(" / ")
   f_env$"^" <- binary_op("^")
   f_env$"[" <- binary_op("_")
   f_env$"==" <- binary_op("=")
   f_env$"<=" <- binary_op(" \\leq ")
   f_env$">=" <- binary_op(" \\geq ")
   f_env$"&" <- binary_op(" \\,\\&\\, ")
   f_env$"|" <- binary_op(" \\mid ")
   f_env$"," <- binary_op(" \\, ")
   # f_env$"," <- function(a, b) {
          # paste0( a, ",", b)
     # }   
   
   # Grouping
   f_env$"{" <- unary_op("\\left{ ", " \\right}")
   f_env$"(" <- unary_op("\\left( ", " \\right)")
   
   # Math functions
   f_env$sin  <- unary_op("\\sin(", ")")
   f_env$cos  <- unary_op("\\cos(", ")")
   f_env$tan  <- unary_op("\\tan(", ")")
   f_env$exp  <- unary_op("\\exp(", ")")
   f_env$expm1  <- unary_op("\\left(\\exp(", ")-1\\right)")
   f_env$sqrt <- unary_op("\\sqrt{", "}")
   f_env$log  <- unary_op("\\log(", ")")
   f_env$log1p  <- unary_op("\\log\\left(1+","\\right)")
   f_env$abs  <- unary_op("\\left| ", "\\right| ")
   f_env$sign  <- unary_op("\\mathop{\\mathrm{sgn}}(", ")")
   f_env$"/"  <- function(a, b) {
          paste0("\\frac{", a, "}{", b, "}")
     }
   f_env$"P"  <- unary_op(" \\mathsf{P}(",")")
   f_env$"F"  <- unary_op(" \\mathsf{F}(",")")
   f_env$"f"  <- unary_op(" \\mathsf{f}(",")")
   f_env$"S"  <- unary_op(" \\mathsf{S}(",")")
   f_env$"H"  <- unary_op(" \\mathsf{H}(",")")
   f_env$"h"  <- unary_op(" \\mathsf{h}(",")")
   f_env$"E"  <- unary_op(" \\mathsf{E}(",")")
   f_env$"V"  <- unary_op(" \\mathsf{V}(",")")  
   f_env$"COV"  <- unary_op(" \\mathsf{COV}(",")")    

   clone_env <- function(env, parent = parent.env(env)) {
       list2env(as.list(env), parent = parent)
    }
   latex_env <- function(expr) {
     greek_env
    }
   # core function
   if (class(object) == "data.frame"){
    cat("%%% LaTeX table generated in R",strsplit(version[['version.string']], ' ')[[1]][3],"by TEX.sde() method",
        "\n")
    cat("%%% Copy and paste the following output in your LaTeX file",
        "\n\n")
    cat(knitr::kable(object, format = "latex",...),
        "\n")
    }else if (class(object) == "MCM.sde"){
    tab <- object$MC
    expr <- parse(text = rownames(tab))
    # greek_test <- as.list(c(mem_sde,mem_sde_tex,expr,
	           # greek,greek1,greek2,greek3,greek4,greek5,
               # greek6,greek7,greek8,greek9,greek10,greek_list,
               # greek_list1,greek_list2,greek_list3,greek_list4,greek_list5,
               # greek_list6,greek_list7,greek_list8,greek_list9,greek_list10))
    names <- all.names(expr )
	greek_env   <- list2env(as.list(c(greek_list,greek_list1,greek_list2,greek_list3,greek_list4,greek_list5,
                                     greek_list6,greek_list7,greek_list8,greek_list9,greek_list10,
                                     mem_sde_list)), parent = emptyenv())
    symbol_list   <- setNames(as.list(names), names)
    symbol_env  <- list2env(symbol_list, parent = f_env)
    greek_env   <- clone_env(greek_env, parent = symbol_env)
    rownames(tab)  <- sapply(1:length(expr),function(i) eval(expr[i], latex_env(expr[i])))
	greek_test <- as.list(c(mem_sde,mem_sde_tex,
	           greek,greek1,greek2,greek3,greek4,greek5,
               greek6,greek7,greek8,greek9,greek10,greek_list,
               greek_list1,greek_list2,greek_list3,greek_list4,greek_list5,
               greek_list6,greek_list7,greek_list8,greek_list9,greek_list10))
    rownames(tab)  <- sapply(1:length(expr),function(i) ifelse(rownames(tab)[i]%in%greek_test,paste0("$", rownames(tab)[i],"$") ,rownames(tab)[i]) )
    colnames(tab)[length(names(tab))] <- "CI( 2.5 \\% , 97.5 \\% )"
    cat("%%% LaTeX table generated in R",strsplit(version[['version.string']], ' ')[[1]][3],"by TEX.sde() method",
        "\n")
    cat("%%% Copy and paste the following output in your LaTeX file",
        "\n\n")
    cat(knitr::kable(tab, format = "latex", escape = FALSE,...),
        "\n")
    }else if (class(object) == "expression"){
    expr <- object
    names <- all.names(expr)
    symbol_list <- setNames(as.list(names), names)
    symbol_env <- list2env(symbol_list, parent = f_env)
    greek_env  <- clone_env(greek_env, parent = symbol_env)
    if (length(expr)== 2){
    dr = eval(expr[[1]], latex_env(expr[[1]]))
    df = eval(expr[[2]], latex_env(expr[[2]]))         
    body <- paste("%%% LaTeX equation generated in R",strsplit(version[['version.string']], ' ')[[1]][3],"by TEX.sde() method")          
    body <- c(body,paste("%%% Copy and paste the following output in your LaTeX file\n"))
    body <- c(body, paste("\\begin{equation}\\label{eq:}"))
    body <- c(body,paste("dX_{t} =",dr,"\\:dt + ",df,"\\:dW_{t}"))
    body <- c(body, paste("\\end{equation}"))
    structure(body, class = "Latex")} else 
    if (length(expr)== 4){
    dr1 = eval(expr[[1]], latex_env(expr[[1]]))
    dr2 = eval(expr[[2]], latex_env(expr[[2]]))
    df1 = eval(expr[[3]], latex_env(expr[[3]]))
    df2 = eval(expr[[4]], latex_env(expr[[4]]))          
    body <- paste("%%% LaTeX equation generated in R",strsplit(version[['version.string']], ' ')[[1]][3],"by TEX.sde() method")          
    body <- c(body,paste("%%% Copy and paste the following output in your LaTeX file\n"))
    body <- c(body, paste("\\begin{equation}\\label{eq:}"))
    body <- c(body, paste("\\begin{cases}"))
    body <- c(body, paste("\\begin{split}"))
    body <- c(body,paste("dX_{t} &=",dr1,"\\:dt + ",df1,"\\:dW_{1,t}","\\\\"))
    body <- c(body,paste("dY_{t} &=",dr2,"\\:dt + ",df2,"\\:dW_{2,t}"))
    body <- c(body, paste("\\end{split}"))
    body <- c(body, paste("\\end{cases}"))
    body <- c(body, paste("\\end{equation}"))
    structure(body, class = "Latex")} else
    if (length(expr)== 6){
    dr1 = eval(expr[[1]], latex_env(expr[[1]]))
    dr2 = eval(expr[[2]], latex_env(expr[[2]]))
    dr3 = eval(expr[[3]], latex_env(expr[[3]]))
    df1 = eval(expr[[4]], latex_env(expr[[4]]))
    df2 = eval(expr[[5]], latex_env(expr[[5]]))
    df3 = eval(expr[[6]], latex_env(expr[[6]]))
    body <- paste("%%% LaTeX equation generated in R",strsplit(version[['version.string']], ' ')[[1]][3],"by TEX.sde() method")          
    body <- c(body,paste("%%% Copy and paste the following output in your LaTeX file\n"))
    body <- c(body, paste("\\begin{equation}\\label{eq:}"))
    body <- c(body, paste("\\begin{cases}"))
    body <- c(body, paste("\\begin{split}"))
    body <- c(body,paste("dX_{t} &=",dr1,"\\:dt + ",df1,"\\:dW_{1,t}","\\\\"))
    body <- c(body,paste("dY_{t} &=",dr2,"\\:dt + ",df2,"\\:dW_{2,t}","\\\\"))
    body <- c(body,paste("dZ_{t} &=",dr3,"\\:dt + ",df3,"\\:dW_{3,t}"))
    body <- c(body, paste("\\end{split}"))
    body <- c(body, paste("\\end{cases}"))
    body <- c(body, paste("\\end{equation}"))
    structure(body, class = "Latex")}
    }else if (class(object) == "MEM.sde"){
    if (object$dim==1){
    expr <- c(parse(text = deparse(object$Means)),parse(text = deparse(object$Var)))
    names <- all.names(expr)
    symbol_list <- setNames(as.list(names), names)
    symbol_env <- list2env(symbol_list, parent = f_env)
    greek_env  <- clone_env(greek_env, parent = symbol_env)
    dm = eval(expr[[1]], latex_env(expr[[1]]))
    dS = eval(expr[[2]], latex_env(expr[[2]]))         
    body <- paste("%%% LaTeX equation generated in R",strsplit(version[['version.string']], ' ')[[1]][3],"by TEX.sde() method")          
    body <- c(body,paste("%%% Copy and paste the following output in your LaTeX file\n"))
    body <- c(body, paste("\\begin{equation}\\label{eq:}"))
    body <- c(body, paste("\\begin{cases}"))
    body <- c(body, paste("\\begin{split}"))
    body <- c(body,paste("\\frac{d}{dt} m(t) &=",dm,"\\\\"))
    body <- c(body,paste("\\frac{d}{dt} S(t) &=",dS))
    body <- c(body, paste("\\end{split}"))
    body <- c(body, paste("\\end{cases}"))
    body <- c(body, paste("\\end{equation}"))
    structure(body, class = "Latex")
   } else if (object$dim==2){
    expr <- c(parse(text = deparse(object$Means[[1]])),parse(text = deparse(object$Means[[2]])),
              parse(text = deparse(object$Var[[1]])),parse(text = deparse(object$Var[[2]])),
              parse(text = deparse(object$Var[[3]])))
   names <- all.names(expr)
   symbol_list <- setNames(as.list(names), names)
   symbol_env <- list2env(symbol_list, parent = f_env)
   greek_env <- clone_env(greek_env, parent = symbol_env)
   dm1 = eval(expr[[1]], latex_env(expr[[1]]))
   dm2 = eval(expr[[2]], latex_env(expr[[2]]))
   dS1 = eval(expr[[3]], latex_env(expr[[3]]))
   dS2 = eval(expr[[4]], latex_env(expr[[4]]))
   dC12 = eval(expr[[5]], latex_env(expr[[5]]))         
   body <- paste("%%% LaTeX equation generated in R",strsplit(version[['version.string']], ' ')[[1]][3],"by TEX.sde() method")          
   body <- c(body,paste("%%% Copy and paste the following output in your LaTeX file\n"))
   body <- c(body, paste("\\begin{equation}\\label{eq:}"))
   body <- c(body, paste("\\begin{cases}"))
   body <- c(body, paste("\\begin{split}"))
   body <- c(body,paste("\\frac{d}{dt} m_{1}(t) ~&=",dm1,"\\\\"))
   body <- c(body,paste("\\frac{d}{dt} m_{2}(t) ~&=",dm2,"\\\\"))
   body <- c(body,paste("\\frac{d}{dt} S_{1}(t) ~&=",dS1,"\\\\"))
   body <- c(body,paste("\\frac{d}{dt} S_{2}(t) ~&=",dS2,"\\\\"))
   body <- c(body,paste("\\frac{d}{dt} C_{12}(t) &=",dC12))
   body <- c(body, paste("\\end{split}"))
   body <- c(body, paste("\\end{cases}"))
   body <- c(body, paste("\\end{equation}"))
   structure(body, class = "Latex")
   } else if (object$dim==3){
    expr <- c(parse(text = deparse(object$Means[[1]])),parse(text = deparse(object$Means[[2]])),parse(text =deparse(object$Means[[3]])),
              parse(text = deparse(object$Var[[1]])),parse(text = deparse(object$Var[[2]])),parse(text = deparse(object$Var[[3]])),
              parse(text = deparse(object$Var[[4]])),parse(text = deparse(object$Var[[5]])),parse(text =deparse(object$Var[[6]])))
   names <- all.names(expr)
   symbol_list <- setNames(as.list(names), names)
   symbol_env <- list2env(symbol_list, parent = f_env)
   greek_env <- clone_env(greek_env, parent = symbol_env)
   dm1 = eval(expr[[1]], latex_env(expr[[1]]))
   dm2 = eval(expr[[2]], latex_env(expr[[2]]))
   dm3 = eval(expr[[3]], latex_env(expr[[3]]))
   dS1 = eval(expr[[4]], latex_env(expr[[4]]))
   dS2 = eval(expr[[5]], latex_env(expr[[5]]))
   dS3 = eval(expr[[6]], latex_env(expr[[6]]))
   dC12 = eval(expr[[7]], latex_env(expr[[7]]))
   dC13 = eval(expr[[8]], latex_env(expr[[8]]))
   dC23 = eval(expr[[9]], latex_env(expr[[9]]))       
   body <- paste("%%% LaTeX equation generated in R",strsplit(version[['version.string']], ' ')[[1]][3],"by TEX.sde() method")          
   body <- c(body,paste("%%% Copy and paste the following output in your LaTeX file\n"))
   body <- c(body, paste("\\begin{equation}\\label{eq:}"))
   body <- c(body, paste("\\begin{cases}"))
   body <- c(body, paste("\\begin{split}"))
   body <- c(body,paste("\\frac{d}{dt} m_{1}(t) ~&=",dm1,"\\\\"))
   body <- c(body,paste("\\frac{d}{dt} m_{2}(t) ~&=",dm2,"\\\\"))
   body <- c(body,paste("\\frac{d}{dt} m_{3}(t) ~&=",dm3,"\\\\"))
   body <- c(body,paste("\\frac{d}{dt} S_{1}(t) ~&=",dS1,"\\\\"))
   body <- c(body,paste("\\frac{d}{dt} S_{2}(t) ~&=",dS2,"\\\\"))
   body <- c(body,paste("\\frac{d}{dt} S_{3}(t) ~&=",dS3,"\\\\"))
   body <- c(body,paste("\\frac{d}{dt} C_{12}(t) &=",dC12,"\\\\"))
   body <- c(body,paste("\\frac{d}{dt} C_{13}(t) &=",dC13,"\\\\"))
   body <- c(body,paste("\\frac{d}{dt} C_{23}(t) &=",dC23))
   body <- c(body, paste("\\end{split}"))
   body <- c(body, paste("\\end{cases}"))
   body <- c(body, paste("\\end{equation}"))
   structure(body, class = "Latex")
   }
   }else {return(paste0("TEX.sde() function not available for this class."))}
}
