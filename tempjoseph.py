

#####
#Supplements to cleanitallup.py
#####

def demoWastedVotes(state, district, demo, party1, party2):
    
    subframe = blockstats.ix[blockstats.ID.isin(list(state.key[state.value == district])), [demo, party1, party2]]
    
    p1Votes  = sum(subframe[party1])
    p2Votes  = sum(subframe[party2])
    numVotes = p1Votes + p2Votes
    
    if p1Votes > p2Votes:
        #p1 wins, waste sum(min/pop*p2)
        #         waste (p1Votes - 0.5*numVotes)*sum(min)/sum(pop)
        return float(sum(blockstats[demo]*blockstats[party2]))/numVotes + \
               (p1Votes - 0.5*numVotes)*sum(blockstats[demo])/numVotes
    else:
        #then p2 wins, do the opposite
        return float(sum(blockstats[demo]*blockstats[party1]))/numVotes + \
               (p2Votes - 0.5*numVotes)*sum(blockstats[demo])/numVotes

def demoEfficiency(state, demo, popcol, party1, party2):
    wasted = [0,0]
    for district in range(ndistricts):
        
        subframe = blockstats.ix[blockstats.ID.isin(list(state.key[state.value == district])), [demo, party1, party2]]
        
        p1Votes  = sum(subframe[party1])
        p2Votes  = sum(subframe[party2])
        numVotes = p1Votes + p2Votes
        
        if p1Votes > p2Votes:
            #p1 wins, waste sum(min/pop*p2)
            #         waste (p1Votes - 0.5*numVotes)*sum(min)/sum(pop)
            wasted[0] = wasted[0] + float(sum(blockstats[demo]*blockstats[party2]))/numVotes + \
                                    (p1Votes - 0.5*numVotes)*sum(blockstats[demo])/numVotes
            wasted[1] = wasted[1] + float(sum((blockstats[popcol] - blockstats[demo])*blockstats[party2]))/numVotes + \
                                    (p1Votes = 0.5*numVotes)*sum((blockstats[popcol] - blockstats[demo]))/numVotes
        else:
            #then p2 wins, do the opposite
            wasted[0] = wasted[0] + float(sum(blockstats[demo]*blockstats[party1]))/numVotes + \
                                    (p2Votes - 0.5*numVotes)*sum(blockstats[demo])/numVotes
            wasted[1] = wasted[1] + float(sum((blockstats[popcol] - blockstats[demo])*blockstats[party1]))/numVotes + \
                                    (p2Votes = 0.5*numVotes)*sum((blockstats[popcol] - blockstats[demo]))/numVotes
    return [float(wasted[0])/sum(blockstats[demo]), float(wasted[1])/sum(blockstats[popcol] - blockstats[demo])]

def singleDemoEfficiency(state, demo, popcol = "population"):
    wasted = [0,0]
    for district in range(ndistricts):
        subframe = blockstats.ix[blockstats.ID.isin(list(state.key[state.value == district])), [demo, popcol]]
        numVotes = sum(subframe[popcol])
        posVotes = sum(subframe[demo])
        negVotes = numVotes - posVotes
        
        if posVotes > negVotes:
            wasted[0] = wasted[0] + posVotes - 0.5*numVotes
            wasted[1] = wasted[1] + negVotes
        else:
            wasted[1] = wasted[1] + negVotes - 0.5*numVotes
            wasted[0] = wasted[0] + posVotes
    return [wasted[0]/sum(blockstats[demo]), wasted[1]/sum(blockstats[popcol])]

def efficiency(state, district):
    #returns difference in percentage of votes wasted.  Negative values benefit R.
    subframe = blockstats.loc[blockstats.ID.isin(list(state.key[state.value == district]))]
    rvotes = sum(subframe['repvotes'])
    dvotes = sum(subframe['demvotes'])
    allvotes = rvotes + dvotes
    
    if rvotes > dvotes:
        wastedR = max(rvotes, dvotes) - 0.5*allvotes
        wastedD = min(rvotes,dvotes)
    else:
        wastedD = max(rvotes, dvotes) - 0.5*allvotes
        wastedR = min(rvotes,dvotes)
    
    return wastedR-wastedD 

def demoEfficiency(state, district, demo1, demo2):
    #returns difference in percentage of ineffective and superfluous votes by demographic.
    subframe = blockstats.loc[blockstats.ID.isin(list(state.key[state.value == district]))]
    demo1Vote = sum(subframe[demo1])
    demo2Vote = sum(subframe[demo2])
    allvotes = demo1Vote + demo2Vote
    
    if demo1Vote > demo2Vote:
        wasted1 = max(demo1Vote, demo2Vote) - 0.5*allvotes
        wasted2 = min(demo1Vote, demo2Vote)
    else:
        wasted2 = max(demo1Vote, demo2Vote) - 0.5*allvotes
        wasted1 = min(demo1Vote, demo2Vote)
    
    return (wasted1, wasted2)


#####
#Supplements to setup_stuff.py
#####

#####
#Supplements to setup.py
#####
