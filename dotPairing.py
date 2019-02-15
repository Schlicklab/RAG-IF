array="(((..((.)).().)))"


#print inputDotBracket(array)

def inputDotBracket(array):
  bracketPositions = []
  results = []
  for i, item in enumerate(array):
    if i == 0 and item == ')':
        pass
        #print("Non sense ! Exit")
        break

    if item == '(':
        bracketPositions.append(i)
    elif item ==')':
        if len(bracketPositions) > 0:
            openingPosition = bracketPositions.pop()
            results.append([openingPosition+1,i+1 ]  )
            #print(openingPosition+1, '-->', i+1)
    else:
        pass
        #print('ERROR: Not a bracket. Word is: %s.' % item)
  return results  


#print array
#print inputDotBracket(array)




