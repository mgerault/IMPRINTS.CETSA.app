#' com_protein_loop
#'
#' Function to know if each element of each element of a list is unique or duplicated; and if so, how many times.
#'
#' @param data A named list containing at least two elements which contains character of numeric elements
#'
#' @return A list containing the common and unique element of the whole list
#'
#' @examples
#' # with numeric
#' x <- list("a" = c(1,2,3,5), "b" = c(55,1,8,99,65,3,4), "c" = c(7,0,1), "d" = c(55,1,65,3,5))
#' com_protein_loop(x)
#'
#' # with character, here these are proteins
#' x <- mineCETSAapp::drug_data
#' com_protein_loop(lapply(x$data, function(l) l$id))
#'
#' @export

com_protein_loop <- function(data){
  if(length(data) <= 1){
    return(data)
  }
  all_common <- Reduce(intersect, data) #get the common element for all the element of the list

  lonelyL <- list()
  rep_L <- list()
  for (i in names(data)){
    lonely <- setdiff(data[[i]], unname(unlist(data[!(names(data) %in% i)])))
    lonelyL[[i]] <- lonely #get the unique element of each element of the list


    repl <- data[[i]][which(is.na(match(data[[i]], c(all_common, lonely))))]

    rep_L[[i]] <- repl  #get the duplicates but not common in all element of the list
  }

  all_common <- list("some" = all_common)
  names(all_common) <- paste(names(data), collapse = " & ")

  all_possi <- all_inter(rep_L)          #call all_inter to know in which groups the duplicates are duplicated
  all_possi <- Reduce(append, all_possi) #get the whole list

  all_possi <- all_possi[!unlist(lapply(all_possi, purrr::is_empty))]  #remove the empty ones

  result <- Reduce(append, list(lonelyL, all_possi, all_common)) #get the whole list (all common, duplicates and unique)
  result <- result[!unlist(lapply(result, purrr::is_empty))]  #remove the empty ones

  return(result)
}



### function nb_comb ###
nb_comb <- function(n){
  l <- list()

  for (i in 2:(n-1)){
    l[[paste0("group_of_", i)]] <- combn(n, i) #get the matrix of possible combinations of groups depending on the length of the list
  }

  return(l)
}

### function all_inter ###
all_inter <- function(data){  #here the data are duplicates but not common in all element of the list
  l <- list()
  a <- nb_comb(length(data)) #get the matrix of possible combinations of groups

  for(i in 1:length(a)){
    for(k in 1:ncol(a[[i]])){
      l[[names(a)[i]]][[paste0(paste(names(data)[a[[i]][,k]], collapse = " & "))]] <- Reduce(intersect, data[a[[i]][,k]])
    } #get the duplicates
  }

  if(length(l) > 1){
    for (k in 1:(length(l) - 1)){
      for(i in 1:(length(l) - k)){
        if(!purrr::is_empty(unname(unlist(l[[(length(l) - k + 1)]])))){
          l[[i]] <- lapply(l[[i]], function(x){
            pr_grp <- which(!is.na(match(x, unname(unlist(l[[(length(l) - k + 1)]])))))
            if(!purrr::is_empty(pr_grp)){
              x <- x[-pr_grp]
            }
            else{
              x <- x
            };
            x
          })
        }
      #if duplicated in three groups, remove the information that is duplicated in two groups
      }
    }
  }

  return(l)
}



