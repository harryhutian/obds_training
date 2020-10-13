ggplot(diamonds)
ggplot(diamonds, aes(x = carat, y = price))
ggplot(diamonds, aes(x = carat, y = price)) +
  geom_point()
ggplot(diamonds, aes(x = carat, y = price, colour = cut)) +
  geom_point()
ggplot(diamonds, aes(x = carat, y = price)) +
  geom_point() +
  geom_smooth()
ggplot(diamonds, aes(x = carat, y = price, colour = cut)) +
  geom_point() +
  geom_smooth()
ggplot(diamonds, aes(x = carat, y = price)) +
  geom_point(aes(colour = cut)) +
  geom_smooth()
