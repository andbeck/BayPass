library(dplyr)
library(magrittr)

AA<-data.frame(ma = runif(10))
AA$mb<-1-AA$ma

BB<-data.frame(ma = runif(10))
BB$mb<-1-BB$ma

AA
BB

df<-rbind(AA,BB)
df$code<-rep(c("A","B"), each = 10)

# rbind(
	# c(AA[1,], BB[1,]),
	# c(AA[2,], BB[2,])
	# # etc
	# )

# YES!
bind_cols(AA,BB)

# Does not work
df %>%
	group_by(code) %>%
		do(vals = bind_cols(data.frame(.)))

# does work, nice trick with []
i=1
out<-bind_cols(split(df, df$code))
out[1:(i+2)!=(i+2)]

# magrittr solution with %$% to pass left to right
out <- df %>% 
  group_by(code) %>% 
  do(each_df = select(., -code)) %$%
  bind_cols(each_df)