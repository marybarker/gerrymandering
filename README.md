A Merry Gander at Gerrymandering
==========

Gerrymandering is the act of redrawing district lines to partition voters in such a way that one party gains a political advantage.  While gerrymandering is sometimes obvious, it can be difficult to quantify.  There is much interest in how to objectively determine whether a district has been gerrymandered, and, if we find a problematic districting plan, what do we do about it?

We have a few goals for this project:
 - Generate many possible district plans that meet requirements.
     - Even populations.
     - Relatively un-weird shapes.
     - Voting Rights Act compliant.
 - Determine a method of evaluating a state's districting plan.
     - Use the many generated plans to identify outlier behavior.
     - Identify patterns to use as a heuristic to spot unfairness.

![NC MH run 250 step](https://github.com/marybarker/gerrymandering/blob/master/visuals/cropSampleNC.gif "North Carolina Congressional District Run")

## Metropolis-Hastings

If we have a state space to navigate as well as a metric for determining whether one state is better than another, we can use the Metropolis-Hastings algorithm to find an optimum.  In our problem, the state space consists of every possible assignment of voting precinct to a congressional district.  To avoid confusion, let's refer to "states" in the "state space" sense as "configurations" instead.

From a given configuration, **A**, there are many "neighboring" configurations.  A neighboring configuration, **B**, will have all precints assigned to the same districts as **A**, except for one.  Any configuration that differs from **A** by only one precinct is a neighbor.  The way we will navigate this space is to randomly choose a neighbor to **A**, and to measure each of their goodness values, *G*(**B**) and *G*(**A**).  If *G*(**B**) is greater than *G*(**A**), we move to **B** so that we can evaluate its neighbors and continue.  If it isn't, then we might go to configuration **B** anyway.  This helps us avoid local maxima.

By moving to **B** with probability ![Prob](https://github.com/marybarker/gerrymandering/blob/master/visuals/CodeCogsEqn.gif "Prob"), we can navigate the space as a Markov Chain, which is satisfying for mathematical purposes.

## Experiment Design

There are several factors which go into the goodness function.

 - Even populations.
 - Relatively un-weird shapes.
 - Voting Rights Act compliant.

We measure the even-ness of populations with two functions.  We consider the maximum population difference among districts, but we also consider the total variation from the mean.

The "un-weirdness" of a shape can be determined in a few different ways.  The official, legal language is that a district should be reasonably, "compact," but that term is vague, and has an alternate meaning in math.  We have used two scores to this end.  For each configuration, we find the ratio of each district's perimeter to the perimeter of a circle with the same area.  We consider both the mean and the maximum of these districts' weirdness.

For Voting Rights Act compliance, the strategy is more nuanced.  As per the 2017 North Carolina ruling, race should not be the only factor in determining district lines.  So in order to ensure that we have majority minority districts, which seems to be a requirement of the VRA, we have to consider larger structures, such as clusters of high-minority precincts, or the gradient from precinct to precinct.

For the experiments we have been running recently, we have been creating a grid with different weights of scores from each of these three factors, and showing how the situation develops over time.

## Current Work

We are working on methods for measuring minority representation so that we can evaluate whether our metrics are effective for meeting the goals of the VRA.

We are also working on parallelizing the algorithm so that we can run Markov Chains on NVIDIA graphics cards.  This would significantly decrease the time necessary to produce a large number of possible districting plans to be a statistical basis.
