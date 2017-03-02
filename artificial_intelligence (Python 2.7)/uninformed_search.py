def generalSearch(problem, fringe, search_type, heuristic):
    """
        
        This is a general function for all type of searchs
        Since all of them shares a similar structure
        """
    visited = []
    path_to_current = []
    
    "Push start state to structure"
    if (search_type=="DFS" or search_type=="BFS"):
        fringe.push((problem.getStartState(), path_to_current))
    elif (search_type=="UCS" or search_type == "A*"):
        fringe.push((problem.getStartState(), path_to_current),None)

    "pop things off and record the path"
    while (fringe.isEmpty() == False):
        operating, path_to_current = fringe.pop()
        skip = False
        if problem.isGoalState(operating):
            "if reach goal state, return path"
            return path_to_current
        
        for node in visited:
            if node == operating:
                skip = True

if skip == False:
    successors = problem.getSuccessors(operating)
        visited.append(operating)
            for successor in successors:
                pos, direction, cost = successor
                latest_path = path_to_current + [direction]
                if (search_type == "DFS" or search_type == "BFS"):
                    fringe.push((pos, latest_path))
                elif (search_type == "UCS"):
                    fringe.push((pos, path_to_current + [direction]),problem.getCostOfActions(latest_path))
                elif (search_type == "A*"):
                    fringe.push((pos, path_to_current + [direction]),problem.getCostOfActions(latest_path)+
                                heuristic(pos,problem))

return None
