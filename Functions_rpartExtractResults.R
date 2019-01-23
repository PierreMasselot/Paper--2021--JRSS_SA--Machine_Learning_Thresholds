#########################################################
#          Functions for extracting thresholds from
#                   rpart results
#
#               Author: Pierre Masselot
#
# Content
# - 
# - 
#########################################################

#-------------------------------------------
#              get.box
#-------------------------------------------

# Extract limits of a subset of predictor space defined by a node in rpart

get.box <- function(object, alert.node = NULL, varnames = names(object$ordered))
# object: rpart object.
# alert.node: node corresponding to the subset. Must be the value in rownames of object$frame. If NULL, take the node with the maximum y mean.
# varnames: names of the variables for which to obtain the subset.

# Value: a list containing three elements
#     - box: a 2 x length(varnames) matrix giving the lower and upper bounds of the subset.
#     - support: the number of observations in the subset
#     - yfun: mean value in the subset
#     - node: the node corresponding to the subset
{
    require(rpart)
    stopifnot(inherits(object,"rpart"))
    tree.data <- eval(object$call[[3]])[,varnames]
    node.box <- apply(tree.data, 2, range)
    node_list <- object$frame
    if (nrow(node_list) == 1) return(list(box = node.box, support = NA, yfun = NA))
    if (is.null(alert.node)){
       alert.node <- row.names(node_list)[which.max(node_list$yval)]
    } else {
       stopifnot(alert.node %in% rownames(node_list))
    }
    alert.path <- path.rpart(object, alert.node, print.it=F)[[1]][-1]
    splits <- strsplit(alert.path, "<|>=")
    sides <- substr(alert.path, regexpr("<|>=",alert.path), regexpr("<|>=",alert.path)+1)
    nsp <- length(alert.path)    
    for (j in 1:nsp){
        split.val <- as.numeric(splits[[j]][2])
        split.var <- grep(splits[[j]][1], varnames)
        node.box[(sides[j] == "< ")+1,split.var] <- split.val
    }
    node.stats <- node_list[alert.node,]
    return(list(box = node.box, support = node.stats$n, yfun = node.stats$yval, node = alert.node))
}

#-------------------------------------------
#              add.partition
#-------------------------------------------

# Add the partition from rpart to an existing plot.
# /!\ works only for bidimensional plots

add.partition <- function(object, vars = names(object$ordered)[1:2], leaf.text = c("all","max","none"), segment.pars = list(), text.pars = list(), digits = 0)
# object: rpart object.
# vars: c(x,y) the two variables of the plot. The default is to take the two first variables.
# leaf.text: if not "none", the mean y value at the leaves is written in the subset. If "all", all leaves are written. If "max", only the leaf with the maximum y mean is written. 
# segment.pars: parameters for the function 'segments' which draws the lines.
# text.pars: parameters for function 'text' if leaf.text != 'none'.
# digits: number of digits for the mean value displayed.
{
    leaf.text <- match.arg(leaf.text)
    ff <- object$frame
    node <- as.numeric(row.names(ff))
    max.node <- node[which.max(ff$yval)]
    leaves <- node[ff$var == "<leaf>"]
    paths <- path.rpart(object,leaves, print.it=F)
    npath <- length(paths)
    for (i in 1:npath){
        ipath <- paths[[i]][-1]
        box <- par("usr")
        ni <- length(ipath)
        splits <- strsplit(ipath, "<|>=")
        sides <- substr(ipath, regexpr("<|>=",ipath), regexpr("<|>=",ipath)+1)
        for (j in 1:ni){
            split.val <- as.numeric(splits[[j]][2])
            split.var <- grep(splits[[j]][1], vars)
            seg.coord <- box
            seg.coord[(1:2)+(split.var-1)*2] <- split.val
            segment.pars <- within(segment.pars,{
               x0 <- seg.coord[1]
               x1 <- seg.coord[2]
               y0 <- seg.coord[3]
               y1 <- seg.coord[4]
            })
            do.call(segments,segment.pars)
            box[(split.var-1)*2+(sides[j] == "< ")+1] <- split.val
        }
        if (leaf.text == "all" || (leaf.text == "max" && names(paths)[i] == as.character(max.node))){
           text.pars$x <- mean(box[1:2])
           text.pars$y <- mean(box[3:4])
           text.pars$labels <- round(ff[as.character(leaves[i]),"yval"], digits = digits) 
           do.call(text, text.pars)
        }
        
    }    
}

#-------------------------------------------
#              colpath_to_node
#-------------------------------------------

# For each node, indicates if it is in the path to the node of interest. 

colpath_to_node <- function(node, tree)
# node: integer, the node of interest. Corresponds to the value in the rownames of tree$frame.
# tree: rpart object.

# Value: a boolean vector where TRUE means that the node is in the path.
{
    path.to.root <- function(node){
        if(node == 1)   # root?
            node
        else            # recurse, %/% 2 gives the parent of node
            c(node, path.to.root(node %/% 2))
    }
    nodes <- as.numeric(row.names(tree$frame))
    res <- nodes %in% path.to.root(as.numeric(node))
    names(res) <- row.names(tree$frame)
    return(res)
}