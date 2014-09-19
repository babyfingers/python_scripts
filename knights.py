knights={"Sir Belvedere":"the Wise", "Sir Lancelot":"the Brave", \
         "Sir Galahad":"the Pure", "Sir Robin":"the Brave", "The Black Knight":"John Clease"}

favorites=knights.keys()
favorites.remove("Sir Robin")
for name, title in knights.iteritems() :
    string = name + ", "
    for fav in favorites :
        if fav == name :
            string += title
            break
    else:
        string += title + ", but not quite so brave as Sir Lancelot."
    print string
