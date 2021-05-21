#' is_in_zone
#'
#' Function to tell if a point is in a delimited zone or not.
#'
#' @param border A dataframe that contains the coordinates of the points that delimit the zone (a square for example)
#' @param target The coordinates of the points you want to check if it is in the zone
#'
#'
#'
#' @return A logical to tell if the point is in the zone
#'
#' @examples
#' square <- data.frame(x = c(1,1,2,2), y = c(1,2,2,1))
#' is_in_zone(square, c(1.5,1.5))
#'
#' @export
#'
is_in_zone <- function(border, target){
  degree = 0
  for (i in 1:nrow(border)) {
    if (i != nrow(border)){
      a = border[i,]
      b = border[i + 1,]
    }
    else{
      a = border[i,]
      b = border[1,]
    }

    # calculate distance of vector
    A = dist(rbind(a,b))[1]
    B = dist(rbind(a,target))[1]
    C = dist(rbind(b, target))[1]

    # calculate direction of vector
    ta_x = a[[1]] - target[1]
    ta_y = a[[2]] - target[2]
    tb_x = b[[1]] - target[1]
    tb_y = b[[2]] - target[2]

    cross = tb_y*ta_x - tb_x*ta_y
    clockwise = cross <= 0

    # calculate sum of angles
    if(clockwise){
      degree = degree + ((acos((B * B + C * C - A * A) / (2.0 * B * C)))*180)/pi
    }
    else{
      degree = degree - ((acos((B * B + C * C - A * A) / (2.0 * B * C)))*180)/pi
    }

  }

  if(abs(round(degree) - 360) <= 3){ #if sum of angles is equal to 360 +/- 3, point is in the zone (3 is an arbitrary confidence parameter)
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}
