
class MinimaxAgent(MultiAgentSearchAgent):

    def getAction(self, gameState):
        possible_actions = []

        # collect all possible actions for current index
        # Since Direction.STOP is always legal, we do not calculate it
        max_action = Directions.STOP
        possible_actions = []

        # collect all possible actions for current index
        # Since Direction.STOP is always legal, we do not calculate it
        for action in gameState.getLegalActions(0):
            if (action != Directions.STOP):
                possible_actions.append(action)

        mini_max_value = float('-inf')
        for action in possible_actions:
            if (mini_max_value < self.miniMax(1, gameState.generateSuccessor(0, action), 1)):
                mini_max_value = self.miniMax(1, gameState.generateSuccessor(0, action), 1)
                max_action = action

        return max_action



    def miniMax(self, depth, state, agentIndex):
        mini_max_value = 0
        if (depth > self.depth) or state.isWin() or state.isLose():
            return self.evaluationFunction(state)
        possible_actions = []

        # collect all possible actions for current index
        # Since Direction.STOP is always legal, we do not calculate it
        for action in state.getLegalActions(agentIndex):
            if (action != Directions.STOP):
                possible_actions.append(action)



        # if it is the last agent, we should go next level to calculate pacman#


        if (agentIndex == state.getNumAgents() - 1):
            next_agent_index = 0
            depth += 1
        else:
            # else calculate next agent
            next_agent_index = agentIndex + 1

        if (agentIndex > 0):
            # if its a ghost, we return the min value
            mini_max_value = float('inf')
            for action in possible_actions:
                if (mini_max_value > self.miniMax(depth, state.generateSuccessor(agentIndex, action), next_agent_index)):
                    mini_max_value = self.miniMax(depth, state.generateSuccessor(agentIndex, action), next_agent_index)

        else:
            # if its the pac man, we return the max value
            mini_max_value = float('-inf')
            for action in possible_actions:
                if (mini_max_value < self.miniMax(depth, state.generateSuccessor(agentIndex, action), next_agent_index)):
                    mini_max_value = self.miniMax(depth, state.generateSuccessor(agentIndex, action), next_agent_index)

        return mini_max_value

class AlphaBetaAgent(MultiAgentSearchAgent):


    def getAction(self, gameState):

        possible_actions = []

        # collect all possible actions for current index
        # Since Direction.STOP is always legal, we do not calculate it
        max_action = Directions.STOP
        possible_actions = []

        # collect all possible actions for current index
        # Since Direction.STOP is always legal, we do not calculate it
        for action in gameState.getLegalActions(0):
            if (action != Directions.STOP):
                possible_actions.append(action)

        alpha = float('-inf')
        beta = float('inf')
        mini_max_value = float('-inf')
        for action in possible_actions:
            if (mini_max_value < self.AB(1, gameState.generateSuccessor(0, action), 1, alpha, beta)):
                mini_max_value = self.AB(1, gameState.generateSuccessor(0, action), 1, alpha, beta)
                max_action = action
                alpha = max(mini_max_value, alpha)

        return max_action



    def AB(self, depth, state, agentIndex, alpha, beta):
        mini_max_value = 0
        if (depth > self.depth) or state.isWin() or state.isLose():
            return self.evaluationFunction(state)
        possible_actions = []

        # collect all possible actions for current index
        # Since Direction.STOP is always legal, we do not calculate i



        # if it is the last agent, we should go next level to calculate pacman#


        if (agentIndex == state.getNumAgents() - 1):
            next_agent_index = 0
            depth += 1
        else:
            # else calculate next agent
            next_agent_index = agentIndex + 1

        if (agentIndex > 0):
            # if its a ghost, we return the min value
            mini_max_value = float('inf')
            for action in state.getLegalActions(agentIndex):
                if action != Directions.STOP:
                    value = self.AB(depth, state.generateSuccessor(agentIndex, action), next_agent_index, alpha, beta)
                    if (mini_max_value > value):
                        mini_max_value = value
                        if (mini_max_value < alpha):
                            return mini_max_value
                        beta = min(mini_max_value, beta)


        else:
            # if its the pacman, we return the max value
            mini_max_value = float('-inf')
            for action in state.getLegalActions(agentIndex):
                if action != Directions.STOP:
                    value = self.AB(depth, state.generateSuccessor(agentIndex, action), next_agent_index, alpha, beta)
                    if (mini_max_value < value):
                        mini_max_value = value
                        if (mini_max_value > beta):
                            return mini_max_value
                        alpha = max(mini_max_value, alpha)

        return mini_max_value





class ExpectimaxAgent(MultiAgentSearchAgent):

    def getAction(self, gameState):
        possible_actions = []

        # collect all possible actions for current index
        # Since Direction.STOP is always legal, we do not calculate it
        max_action = Directions.STOP
        possible_actions = []

        # collect all possible actions for current index
        # Since Direction.STOP is always legal, we do not calculate it
        for action in gameState.getLegalActions(0):
            if (action != Directions.STOP):
                possible_actions.append(action)

        mini_max_value = float('-inf')
        for action in possible_actions:
            if (mini_max_value < self.miniMax(1, gameState.generateSuccessor(0, action), 1)):
                mini_max_value = self.miniMax(1, gameState.generateSuccessor(0, action), 1)
                max_action = action

        return max_action



    def miniMax(self, depth, state, agentIndex):
        mini_max_value = 0
        if (depth > self.depth) or state.isWin() or state.isLose():
            return self.evaluationFunction(state)
        possible_actions = []

        # collect all possible actions for current index
        # Since Direction.STOP is always legal, we do not calculate it
        for action in state.getLegalActions(agentIndex):
            if (action != Directions.STOP):
                possible_actions.append(action)



        # if it is the last agent, we should go next level to calculate pacman#


        if (agentIndex == state.getNumAgents() - 1):
            next_agent_index = 0
            depth += 1
        else:
            # else calculate next agent
            next_agent_index = agentIndex + 1

        if (agentIndex > 0):
            # if its a ghost, we return the min value
            total = 0
            for action in possible_actions:
                total += self.miniMax(depth, state.generateSuccessor(agentIndex, action), next_agent_index)

            mini_max_value = total/len(possible_actions)

        else:
            # if its the pac man, we return the max value
            mini_max_value = float('-inf')
            for action in possible_actions:
                if (mini_max_value < self.miniMax(depth, state.generateSuccessor(agentIndex, action), next_agent_index)):
                    mini_max_value = self.miniMax(depth, state.generateSuccessor(agentIndex, action), next_agent_index)

        return mini_max_value

def betterEvaluationFunction(currentGameState):

    # In this evaluation, I first calculate the manhattanDistance to the nearest
    # food, then I find the distance to the closed ghost. Finally I write
    # an algorithm to let pacman only eat ghost when it has enough confidence, else
    # the pacman will avoid ghosts

    food_pos = currentGameState.getFood()
    pac_pos = currentGameState.getPacmanPosition()
    ghost_state = currentGameState.getGhostStates()
    currFood = currentGameState.getFood()
    mininum_food_dist = -1
    # find the minimum dist to next food
    for food in currFood.asList():
        if (mininum_food_dist == -1):
            mininum_food_dist = manhattanDistance(food, pac_pos)
        elif (manhattanDistance(food, pac_pos) < mininum_food_dist):
            mininum_food_dist = manhattanDistance(food, pac_pos)

    mininum_ghost_dist = 999999
    if (len(ghost_state) > 0):
        nearest_ghost_state = ghost_state[0]
    # find the minimum dist to nearest ghost for currentState
    for ghost_s in ghost_state:
        if (mininum_ghost_dist == 999999):
            mininum_ghost_dist = manhattanDistance(ghost_s.getPosition(), pac_pos)
            nearest_ghost_state = ghost_s
        elif (manhattanDistance(ghost_s.getPosition(), pac_pos) < mininum_ghost_dist):
            mininum_ghost_dist = manhattanDistance(ghost_s.getPosition(), pac_pos)
            nearest_ghost_state = ghost_s

    bonus = 0

    if (mininum_ghost_dist < 3):
        if (nearest_ghost_state.scaredTimer > 3):
            bonus = 5
        else:
            bonus = -5

    return currentGameState.getScore() - mininum_food_dist + bonus


# Abbreviation
better = betterEvaluationFunction

class ContestAgent(MultiAgentSearchAgent):

    def getAction(self, gameState):
        util.raiseNotDefined()

