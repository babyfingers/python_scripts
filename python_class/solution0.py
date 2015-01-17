elements = {1 : "H", 2 : "He", 3 : "Li", 4 : "Be", 5: "B", 6 : "C", 7: "N", 8 :"O", 9 : "F", 10 : "Ne"}
userInput = ""
while (userInput != 0):
    userInput = input("Please enter an atomic number.\n")
    if userInput in elements:
        print "Your input corresponds to " + elements[userInput] +"."
    elif (userInput == 0):
        break
    else:
        print "Sorry, I don't understand."


