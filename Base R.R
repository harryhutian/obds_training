zoo <- c("giraffe","tiger","monkey","zebra","lion")
attr(zoo, "names") <- c("Ben","George","Harry","Jim","Ken")
names(zoo)
zoo
zoo[c("Ben","Harry")]
length(zoo)
zoo[length(zoo)]
zoo[4] <- c("bear")
zoo
num <- c(6, 3, 2, 18)
bool <- c(TRUE, FALSE, TRUE, TRUE)
sum(num[bool])
num*bool
num_2 <- c(3, 4, 12)
num + num_2
missing <- c(9, 4, NA, 2, 43, NA, 4, 2)
is.na(missing)
sum(is.na(missing))
which(is.na(missing))
signal <- factor(c("red","yellow","yellow","green","red","green"),
levels = c("red","yellow","green"), ordered = TRUE)
levels(signal)
signal[2] <- c("blue")
new_list <- list(zoo, bool, signal, num)
new_list
new_list[2:3]
subset_list <- (new_list[2:3])
class(subset_list)
new_list[c(2,4)]
new_list[[3]][[2]]
class(new_list[[3]])
attr(new_list, "names") <- c("Mon","Tue","Wed","Thu")
new_list
new_list$Wed
class(new_list$Mon)
levels(new_list$Wed)
levels(new_list$Wed) <- c("yellow","red","green")
