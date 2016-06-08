def backtracking(choice_points, choices, assignable):
  N = len(choices)
  M = len(choice_points)
  solutions = {}
  cp = 0
  c = 0
  backtrack = False
  end = False
  
  while(not end):
    # going forward
    while(not backtrack):
      if(assignable(cp, c, solutions)):
        solutions[cp] = c
        if(cp==M-1):
          yield {choice_points[k]:choices[v] for k,v in solutions.iteritems()}
          del solutions[cp]
          if not c==N-1:
            c+=1
          else:
            backtrack = True
        else:
          cp+=1
          c=0
      elif(c!=N-1):
        c+=1
      else:
        backtrack = True
    # going backward
    end = (cp==0)
    while(backtrack and not end):
      cp-=1
      c = solutions.pop(cp)
      if(not c==N-1):
        c+=1
        backtrack = False
      elif(cp==0):
        end = True


    
