#use strict;
use Carp;
use multiTree;
use multiGall;
use multiMethods;

my $matrices_ref = matrix_read_file($ARGV[0]);
my $perfectPhylogenyCount = 0; my $networkRecombCount = 0; my $multiCrossCount = 0;
my $num_matrices = @{$matrices_ref};

my $bMinTotalNumCrossOvers = 1;
if (defined $ARGV[1]) {
	$bMinTotalNumCrossOvers = $ARGV[1];
}

#loop from here...
for (my $j = 0; $j < $num_matrices; $j++) {
	my $matrix_name = $matrices_ref->[$j];
	my $matrix_ref = \@{$matrix_name};

	my $uniqueRow_ref = multiMethods->getUniqueRows($matrix_ref);
	$matrix_ref = multiMethods->getUniqueRowsMatrix($matrix_ref, $uniqueRow_ref);
	my $num_rows = @{$matrix_ref};

	print "\n\nProcessing matrix \"$matrix_name\"...\n\n";

	print "Unique rows matrix is\n";
	# printing each row of the matrix
	for (my $i = 0; $i < $num_rows; $i++) {
		print @{$matrix_ref->[$i]}, " (row ", $i+1, ")", "\n";
	}
	print "\n";

	 my $bNoConflicts; my $confColHash_ref;
	($bNoConflicts, $confColHash_ref) = findConflicts($matrix_ref);

	my $num_cols =  @{$matrix_ref->[0]};

	#for (my $i = 0; $i < $num_cols; $i++) {
	#	print "Column $i conflicts with the columns @{$confColHash_ref->{$i}}\n";
	#}
	#print "\n";


	if ($bNoConflicts)	{
		print "There are no conflicting columns, so perfect phylogeny exists.\n";
		$perfectPhylogenyCount++;
		$networkRecombCount++;
		next; # skip rest of the for() loop's block
	}

	#print "For the set of species, perfect phylogeny without recombination doesn't exist\n";
	#print "Testing for phylogeny with recombination begins ....\n\n";

	my $ccs_ref = getConnectedComponents($matrix_ref, $confColHash_ref);
	my $num_cc = @{$ccs_ref};

	#print "\nThe $num_cc connected components are as follows.\n";
	#for (my $i = 0; $i < $num_cc; $i++) {
	#	print "@{$ccs_ref->[$i]}\n";
	#}
	#print "\n";


	# MULTIPLE CROSS-OVER. NO SPECIAL TREATMENT FOR DUPLICATE COLUMNS

	my $ccMatrices_ref = getCCMatrices($matrix_ref, $ccs_ref);

	my $splitMatrix_ref; my $splitMatrixColumnMap_ref;
	($splitMatrix_ref, $splitMatrixColumnMap_ref) = getSplitMatrix($ccMatrices_ref);

	#print "The final matrix of splits is\n";
	#for (my $j = 0; $j < $num_rows; $j++) {
	#	print @{$splitMatrix_ref->[$j]}, " (row ", $j+1, ")", "\n";
	#}
	#print "\n";

	#print "The split columns map to\n";
	#foreach my $item ( @{$splitMatrixColumnMap_ref} ) {
	#	print "$item->[0], @{$item->[1]}\n";
	#}
	#print "\n";

	# relabel original matrix with the first leaf as root
	$splitMatrix_ref = relabelMatrix($splitMatrix_ref, 0);

	#print "The relabelled matrix of splits is\n";
	#for (my $j = 0; $j < $num_rows; $j++) {
	#	print @{$splitMatrix_ref->[$j]}, " (row ", $j+1, ")", "\n";
	#}
	#print "\n";

	# construct T_bar from the matrix of splits
	my $rootNode = multiMethods->buildPerfPhylogeny($splitMatrix_ref);
	if ($rootNode == -1) {
		print "Matrix of splits not forming perfect phylogeny; incorrect algorithm?\n";
		next;
		#croak "Matrix of splits not forming perfect phylogeny; incorrect algorithm?\n";
	}

	#print "T_bar ...\n";
	#$rootNode->showSubtree();
	#print "\n";
	$rootNode->reLabelMultipleSites($splitMatrixColumnMap_ref); # remap T_bar edge labels to splits
	#$rootNode->showSubtreeCompoundEdgeLabels();
	#print "\n";

	# get a hash of nodes to be expanded into galls for all non-trivial CCs
	my %gallNodesHash;
	$rootNode->getGallNodes(\%gallNodesHash, $num_cc);

	# check if all non-trivial ccs have an expandable node and that trivial ccs don't have any in T_bar
	for (my $i = 0; $i < $num_cc; $i++) {
		if ( scalar(@{$ccs_ref->[$i]}) > 1 ) {
			if ( not defined $gallNodesHash{$i} ) {
				croak "Non-trivial CC @{$ccs_ref->[$i]} does not have an expandable node in T_bar. Something wrong!\n";
			}
		}
		else {
			if (defined $gallNodesHash{$i} ) {
				croak "Trivial CC @{$ccs_ref->[$i]} has an expandable node in T_bar. Something wrong!\n";
			}
		}
	}

	#print "Nodes to become galls are...\n";
	#foreach my $i (keys % gallNodesHash) {
	#	my $node = $gallNodesHash{$i};
	#	print "CC# $i: ";
  	#	$node->showNodeCompoundEdgeLabels();
	#}

	my $noGallCC = 0;
	my %ccGallableTreesHash;

	# for each non-trivial cc node, check if expandable into gall and direct edges in T_bar accordingly
	foreach my $ccNum (keys % gallNodesHash) {
		my $cc_ref = $ccs_ref->[$ccNum];

		#print "\nFinding gall for non-trivial cc $ccNum...\n";
		my $gallableTrees_ref = getGallStructures($ccMatrices_ref->[$ccNum], $cc_ref);
		if ( scalar(@{$gallableTrees_ref}) == 0 ) {
			print "Non-trivial cc# $ccNum cannot be galled. No galled-tree exists for matrix.\n";
			$noGallCC = 1;
			last;
		}

		$ccGallableTreesHash{$ccNum} = $gallableTrees_ref;
	}

	if ( $noGallCC == 1 ) {
		print "At least one non-trivial cc cannot be galled. No galled-tree solution.\n";
		next;
	}

	#print "\n";

	# T_bar nodes with an immediately branching leaf node can absorb that leaf into themselves
	$rootNode->mergeSelfLeaf();
	#print "After merging the branching leaf nodes...\n";
	#$rootNode->showSubtreeCompoundEdgeLabels();

	# converting compound edge labels of trivial ccs into simple edge labels (site number)
	$rootNode->simplifyTrivialCCEdgeLabels($ccs_ref);
	#print "After simplifying trivial-cc edge labels...\n";
	#$rootNode->showSubtreeCompoundEdgeLabels();

	# direct edges in T_bar when there is only one possible arrangement for a gall
	directRecombEdges(\%gallNodesHash, \%ccGallableTreesHash);

	my @rootableNodes;
	pickRoots($rootNode, \@rootableNodes); # select all possible roots in partially directed tree

	if ( scalar(@rootableNodes) == 0 ) {  # T_bar is not rootable
		print "Could not find any root node in T_bar\n";
		next;
	}

	# Compute the following for every possible root
	# a) Total number of crossovers over all galls
	# b) The max. number of crossovers in a gall over all galls

	my @crossOvers; # stores every possible root and corresponding array with values " [ a), b) ] " as described above

	foreach my $root (@rootableNodes) {
		$root->directEdgesAway(); # Direct all edges away from this possible root
		my @crossOverData;
		getCrossOverData(\%gallNodesHash, \%ccGallableTreesHash, \@crossOverData);
		push @crossOvers, [$root, \@crossOverData];
	}

	if ($bMinTotalNumCrossOvers == 1) {
		# sort the crossOvers by the total number of crossovers over all galls
		@crossOvers = sort {$a->[1]->[0] <=> $b->[1]->[0]} @crossOvers;
	}
	else {
		# sort the crossOvers by the max. number of crossovers in a gall over all galls
		@crossOvers = sort {$a->[1]->[1] <=> $b->[1]->[1]} @crossOvers;
	}

	$rootNode = $crossOvers[0]->[0]; # the root that minimizes the total number of crossovers or
					 # minimizes the max. number of crossovers in a gall over all galls

 	$rootNode->directEdgesAway(); # newly selected root. Direct all edges away from it
	#print "Newly rooted T_bar is...\n";
	#$rootNode->showSubtreeCompoundEdgeLabels();
	#print "\n";

	croak "Root sequence must be leaf. Bug!\n" if (not defined $rootNode->sequence());
	my $rootSeqID = $rootNode->sequence();

	# expand gall nodes into galls
	expandGalls(\%gallNodesHash, \%ccGallableTreesHash);
	#print "After expanding gall nodes...\n";
	#$rootNode->showSubtreeCompoundEdgeLabels();

	# check if a site occurs at most once in the galled tree corresponding to unique columns matrix (which would
	# imply that it is also true for the galled tree corresponding to the original matrix with duplicate columns)
	my @siteCountVector = (0..$num_cols-1);
	for (my $i = 0; $i < scalar(@siteCountVector); $i++) {
		$siteCountVector[$i] = 0; # initialize the counters to 0 for each site
	}

	# check if every leaf occurs exactly once in the galled tree corresponding to unique rows matrix
	my @leafCountVector = (0..$num_rows-1);
	for (my $i = 0; $i < scalar(@leafCountVector); $i++) {
		$leafCountVector[$i] = 0; # initialize the counters to 0 for each leaf
	}

	my $onSites_ref = [];

	# used to find the zero-One state (like 0101100011) at each node. Can't just refer back to
	# sequences in matrix since some nodes have no corresponding sequence in the matrix.
	# Also checks if the sequences generated match the input sequences. returns -1 if failed.
	if ( $rootNode->findOnSites($num_cols, $onSites_ref, $matrix_ref,
					$rootSeqID, \@siteCountVector, \@leafCountVector) == -1 ) {
		print "Computed sequences and alleged corresponding input sequences are not equal\n";
		next;
	}

	#print "leaf count vector = @leafCountVector\n";
	if ( multiMethods->leafExactlyOnce(\@leafCountVector) == 0 ) {
		print "One or more leaves missing or occuring more than once. Not a galled tree.\n";
		next;
	}

	#print "Site count vector = @siteCountVector\n";
	if ( multiMethods->siteAtmostOnce(\@siteCountVector) == 0 ) {
		print "One or more sites mutating more than once. Not a galled tree.\n";
		next;
	}

	# at this stage, the network has been found
	$networkRecombCount++;

	$multiCrossCount++ if ( $rootNode->hasMultipleCrossovers() == 1 );


	print "The galled tree is...\n";
	$rootNode->showSubtreeCompoundEdgeLabels();
}


$perfandsingle = $networkRecombCount - $multiCrossCount;
print "\n\nSummary: $networkRecombCount/$num_matrices of the input matrices can be derived on a galled-tree.\n";
print "$multiCrossCount of these galled-trees require(s) at least one multiple-crossover recombination,
and the other $perfandsingle of the solutions are either 
perfect phylogenies or require only single crossover recombination.\n";

# expands a gallable node in T_bar into a gall
sub expandGallNode {
	return unless (@_ == 6);

	my $gallNode = $_[0];
	my $gallableTree_ref = $_[1];
	my $coalescentSeq_ref = $_[2];
	my $recombSeq_ref = $_[3];
	my $ccNum = $_[4];
	my $crossOverPts_ref = $_[5];

	# create a new gall
	my $gall = multiGall->new();
	$gall->crossOverPts($crossOverPts_ref);

	# find the coalescent node in the gallable tree
	my $coalescNode = $gallableTree_ref;
	while ( $coalescNode && (isEqual($coalescNode->zeroOneState(), $coalescentSeq_ref) == 0)  ) {
		#print "@{$coalescNode->zeroOneState()}\n";
		$coalescNode = $coalescNode->children()->[0];
	}

	croak "Coalescent node must be present in tree. Bug!\n" if (not defined $coalescNode);

	# create suffix nodes of gall
	my $sufNode = $coalescNode->children()->[0];
	while ($sufNode) {
		my $node = multiTree::Node->new($sufNode->allsites()->[0]);
		$node->allsites($sufNode->allsites()); # simple edge labels
		$node->zeroOneState($sufNode->zeroOneState());
		push @{$gall->suffix()}, $node;

		$sufNode = $sufNode->children()->[0];
	}

	# create prefix nodes of gall
	my $preNode = $coalescNode->parent();
	while ($preNode) {
		my $child = $preNode->children()->[0];
		my $node = multiTree::Node->new($child->allsites()->[0]);
		$node->allsites($child->allsites()); # simple edge labels
		$node->zeroOneState($preNode->zeroOneState());
		push @{$gall->prefix()}, $node;

		$preNode = $preNode->parent();
	}

	croak "Gallable node is interior. Must have parent. Bug!\n" if (not defined $gallNode->parent());

	# process coalescent node
	my $allSites_ref = $gallNode->allsites(); # compound edge labels
	if ( scalar(@{$allSites_ref}) == 1 ) { # no site mutations just above the coalescent
		 # parent of gall node must point to gall instead of gall node
		$gallNode->parent()->addChild($gall); # updates gall's parent also
	}
	else { # sites mutating above coalescent. create a new coalescent node
		for (my $i = 0; $i < scalar(@{$allSites_ref}); $i++) {
			my $label = $allSites_ref->[$i];
			my $gallCCnum = $label->[0];
			if ($ccNum == $gallCCnum) {
				splice @{$allSites_ref}, $i, 1;
				$gallNode->site($allSites_ref->[0]);
				last;
			}
		}

		my $newNode = multiTree::Node->new($allSites_ref->[0]); # compound site label
		$newNode->allsites($allSites_ref); # no sequence for newly created coalescent node
		$gallNode->parent()->addChild($newNode); # updates new node's parent also
		$newNode->addChild($gall);
	}

	# gall node is replaced by gall as child of the coalescent; so its no longer a child
	my $children = $gallNode->parent()->children();
	for (my $i = 0; $i < scalar(@{$children}); $i++) {
		my $child = $children->[$i];
		if ($child == $gallNode) {
			splice @{$children}, $i, 1;
			last;
		}
	}

	# process recombination node
	my $recombNode = getGallNodeChild($gallNode, $ccNum, $recombSeq_ref);
	croak "Recombination node must exist in T_bar. Bug!\n" if ($recombNode == 0);

	$gall->zeroOneState($recombSeq_ref);
	$allSites_ref = $recombNode->allsites();
	if ( scalar(@{$allSites_ref}) == 1 ) { # only immediate child below recombination node
		$gall->sequence($recombNode->sequence());
		my $recombChildren = $recombNode->children();
		for (my $i = 0; $i < scalar(@{$recombChildren}); $i++) {
			$gall->addChild($recombChildren->[$i]); # ensures parent->child and child->parent links
		}
	}
	else { # no immediate child. recombination node undergoes mutation(s) before leading to leaves
		# gall's sequence is undefined
		for (my $i = 0; $i < scalar(@{$allSites_ref}); $i++) {
			my $label = $allSites_ref->[$i];
			my $gallCCnum = $label->[0];
			if ($ccNum == $gallCCnum) {
				splice @{$allSites_ref}, $i, 1;
				$recombNode->site($allSites_ref->[0]);
				last;
			}
		}
		$gall->addChild($recombNode);
	}

	# process prefix/suffix nodes
	for (my $side = 0; $side < 2; $side++) {
		my $nodes_ref = ($side == 0) ? $gall->prefix() : $gall->suffix();
		my $numNodes = @{$nodes_ref};

		for (my $j = 0; $j < $numNodes; $j++) {
			my $psNode = $nodes_ref->[$j]; # prefix/suffix node

			my $node = getGallNodeChild($gallNode, $ccNum, $psNode->zeroOneState());
			croak "Prefix/Suffix node must exist in T_bar. Bug!\n" if ($node == 0);

			$allSites_ref = $node->allsites();
			if ( scalar(@{$allSites_ref}) == 1 ) { # only immediate child below prefix/suffix node
				$psNode->sequence($node->sequence());
				my $psChildren = $node->children();
				for (my $i = 0; $i < scalar(@{$psChildren}); $i++) {
					$psNode->addChild($psChildren->[$i]); # ensures parent->child and child->parent links
				}
			}
			else { # no immediate child. prefix/suffix node undergoes mutation(s) before leading to leaves
				# prefix/suffix node's sequence is undefined
				for (my $i = 0; $i < scalar(@{$allSites_ref}); $i++) {
					my $label = $allSites_ref->[$i];
					my $gallCCnum = $label->[0];
					if ($ccNum == $gallCCnum) {
						splice @{$allSites_ref}, $i, 1;
						$node->site($allSites_ref->[0]);
						last;
					}
				}
				$psNode->addChild($node);
			}
		}
	}


}

# returns child node of gallable node in T_bar with the cc restricted sequence
sub getGallNodeChild {
	return unless (@_ == 3);

	my $gallNode = $_[0];
	my $ccNum = $_[1];
	my $ccRestrictSeq = $_[2];

	my $childNode = 0;

	my $children = $gallNode->children();
	for (my $i = 0; $i < scalar(@{$children}); $i++) {
		my $child = $children->[$i];
		my $allSites_ref = $child->allsites();
		for (my $i = 0; $i < scalar(@{$allSites_ref}); $i++) {
			my $label = $allSites_ref->[$i];
			my $gallCCnum = $label->[0];
			my $ccSeq = $label->[1];
			if ( ($ccNum == $gallCCnum) && (isEqual($ccSeq, $ccRestrictSeq) == 1) ) {
				$childNode = $child;
				last;
			}
		}
		last if ($childNode != 0);
	}

	return $childNode;
}


# compute the total number of crossovers in galled-tree, max. # of crossovers in a gall over all the galls
sub getCrossOverData {
	return unless (@_ == 3);

	my $gallNodesHash_ref = $_[0];
	my $ccGallableTreesHash_ref = $_[1];
	my $crossOverData_ref = $_[2];

	my @numCrossOvers;
	my $totalNumCrossOvers = 0;
	my $maxCrossOversGall = -1; # initializing to some convenient value

	foreach my $ccNum (keys % {$ccGallableTreesHash_ref}) {
		my $ccGallNode_ref = $gallNodesHash_ref->{$ccNum};
		my $coalescentSeq_ref = 0;
		my @allsites = @{$ccGallNode_ref->allsites()};

		croak "Gallable node must have incoming edge with cc label. Implementation bug!\n" if ( scalar(@allsites) == 0 );

		for (my $i = 0; $i < scalar(@allsites); $i++) {
			my $item = $allsites[$i];
			if ($item->[0] == $ccNum) { # edge label corresponding to the non-trivial cc
				$coalescentSeq_ref = $item->[1]; # incoming node defines the gall's coalescent node
				last;
			}
		}

		croak "Missing coalescent sequence. Implementation bug!\n" if ($coalescentSeq_ref == 0);
		#print "Coalescent sequence for cc# ", $ccNum, " is @{$coalescentSeq_ref}\n";

		my $gallableTrees_ref = $ccGallableTreesHash_ref->{$ccNum};

		my @numCrossOversCC; # collect all possible # of crossovers for this gall (over all solutions)

		for (my $i = 0; $i < scalar(@{$gallableTrees_ref}); $i++) {
			my $gallableSolution = $gallableTrees_ref->[$i];
			my $deletedSeq_ref = $gallableSolution->[1];
			my $crossOverPts_ref = $gallableSolution->[2];

			if (isEqual($coalescentSeq_ref, $deletedSeq_ref) == 0) { # coalescent node is not the recombination node
				push @numCrossOversCC, scalar(@{$crossOverPts_ref});
			}
		}

		croak "Missing gallable tree. Implementation bug!\n" if (scalar(@numCrossOversCC) == 0);

		my $minNumCrossOversCC = $numCrossOversCC[0]; # compute the min # of crossovers for this gall

		for (my $i = 1; $i < scalar(@numCrossOversCC); $i++) {
			$minNumCrossOversCC = $numCrossOversCC[$i] if ($numCrossOversCC[$i] < $minNumCrossOversCC);
		}

		push @numCrossOvers, $minNumCrossOversCC;
		$totalNumCrossOvers += $minNumCrossOversCC;
		$maxCrossOversGall = $minNumCrossOversCC if ($minNumCrossOversCC > $maxCrossOversGall);
	}

 	@{$crossOverData_ref} = ($totalNumCrossOvers, $maxCrossOversGall);

}


# expands all gallable nodes in T_bar into their corresponding galls
sub expandGalls {
	return unless (@_ == 2);

	my $gallNodesHash_ref = $_[0];
	my $ccGallableTreesHash_ref = $_[1];

	foreach my $ccNum (keys % {$ccGallableTreesHash_ref}) {
		my $ccGallNode_ref = $gallNodesHash_ref->{$ccNum};
		my $coalescentSeq_ref = 0;
		my @allsites = @{$ccGallNode_ref->allsites()};

		croak "Gallable node must have incoming edge with cc label. Implementation bug!\n" if ( scalar(@allsites) == 0 );

		for (my $i = 0; $i < scalar(@allsites); $i++) {
			my $item = $allsites[$i];
			if ($item->[0] == $ccNum) { # edge label corresponding to the non-trivial cc
				$coalescentSeq_ref = $item->[1]; # incoming node defines the gall's coalescent node
				last;
			}
		}

		croak "Missing coalescent sequence. Implementation bug!\n" if ($coalescentSeq_ref == 0);
		#print "Coalescent sequence for cc# ", $ccNum, " is @{$coalescentSeq_ref}\n";

		my $gallableTrees_ref = $ccGallableTreesHash_ref->{$ccNum};
		my $recombSeq_ref;
		my $gallableTree;

		my @gallSolutions;

		for (my $i = 0; $i < scalar(@{$gallableTrees_ref}); $i++) {
			my $gallableSolution = $gallableTrees_ref->[$i];
			my $tree = $gallableSolution->[0];
			my $deletedSeq_ref = $gallableSolution->[1];
			my $crossOverPts_ref = $gallableSolution->[2];

			if (isEqual($coalescentSeq_ref, $deletedSeq_ref) == 0) { # coalescent node is not the recombination node
				$gallableTree = $tree;
				$recombSeq_ref = $deletedSeq_ref;

				push @gallSolutions, [scalar(@{$crossOverPts_ref}), $gallableTree, $recombSeq_ref, $crossOverPts_ref];
			}
		}

		croak "Missing gallable tree. Implementation bug!\n" if (scalar(@gallSolutions) == 0);

		@gallSolutions = sort {$b->[0] <=> $a->[0]} @gallSolutions; # sort by the number of crossovers

		# solution with the minimum number of crossovers
		$gallableTree = $gallSolutions[0]->[1];
		$recombSeq_ref = $gallSolutions[0]->[2];
		my $crossOverPts_ref = $gallSolutions[0]->[3];

		#print "Gallable tree is ...\n";
		#$gallableTree->showSubtree();

		# replace gallable node with gall
		expandGallNode($ccGallNode_ref, $gallableTree, $coalescentSeq_ref, $recombSeq_ref, $ccNum, $crossOverPts_ref);
	}

}


sub pickRoots {
	return unless (@_ == 2);

	my $root = $_[0];
	my $rootableNodes_ref = $_[1];

	my @allNodes;
	$root->collectAllNodes(\@allNodes);
	my $numNodes = scalar(@allNodes);

	my $nodeNum;
	for (my $nodeNum = 0; $nodeNum < $numNodes; $nodeNum++) {
		my $node_ref = $allNodes[$nodeNum];
		if ($node_ref->isRootable() == 1) {
			#print "A possible root for T_bar is...";
			#$node_ref->showNodeCompoundEdgeLabels();
			#print "\n\n";

			push @{$rootableNodes_ref}, $node_ref;
		}
	}
}


sub directEdge  {
	return unless (@_ == 3);

	my $node_ref = $_[0];
	my $ccNum = $_[1];
	my $ccSeq_ref = $_[2];

	my @parentSites = @{$node_ref->allsites()};
	my $numParentLabels = scalar(@parentSites);


	for (my $i = 0; $i < $numParentLabels; $i++) {
		my $label = $parentSites[$i];
		my $ccLabel = $label->[0];
		my $seqLabel_ref = $label->[1];

		if ( $ccNum == $ccLabel ) {
			if ( isEqual($ccSeq_ref, $seqLabel_ref) == 1 ) { # found edge with matching label = (cc#, ccRestrictedSequence)
				my $parentNode_ref = $node_ref->parent();
				my $parentChildSeekHash_ref = $parentNode_ref->seekChildren();
				croak "Parent not remembering directionality towards child. Impl. bug!" if ( not defined $parentChildSeekHash_ref->{$node_ref} );
				$parentChildSeekHash_ref->{$node_ref} = 0; # parent cannot seek this node
				#print "Directed edge to parent\n";
				return;
			}
		}
	}


	my @childNodes = @{$node_ref->children()};
	my $numChildren = scalar(@childNodes);
	for (my $i = 0; $i < $numChildren; $i++) {
		my $childNode_ref = $childNodes[$i];
		my @childSites = @{$childNode_ref->allsites()};
		my $numChildLabels = scalar(@childSites);

		for (my $i = 0; $i < $numChildLabels; $i++) {
			my $label = $childSites[$i];
			my $ccLabel = $label->[0];
			my $seqLabel_ref = $label->[1];

			if ( $ccNum == $ccLabel ) {
				if ( isEqual($ccSeq_ref, $seqLabel_ref) == 1 ) { # found edge with matching label = (cc#, ccRestrictedSequence)
					$childNode_ref->seekParent(0); # this child cannot seek its parent node
					#print "Directed edge to child\n";
					return;
				}
			}
		}
	}

	croak "Cannot direct edge around gallable node. Implementation bug!\n";
}


sub directRecombEdges {
	return unless (@_ == 2);

	my $gallNodesHash_ref = $_[0];
	my $ccGallableTreesHash_ref = $_[1];

	foreach my $ccNum (keys % {$ccGallableTreesHash_ref}) {
		my $gallableTrees_ref = $ccGallableTreesHash_ref->{$ccNum};

		#print "Collected following trees for non-trivial cc# $ccNum...\n";
		#for (my $i = 0; $i < scalar(@{$gallableTrees_ref}); $i++) {
		#	my $gallableSolution = $gallableTrees_ref->[$i];
		#	my $tree = $gallableSolution->[0];
		#	my $deletedSeq_ref = $gallableSolution->[1];

		#	print "For deleted sequence ", @{$deletedSeq_ref}, " a gallable tree is\n";
		#	$tree->showSubtree();
		#}
		#print "\n";

		# direct edges in T_bar if only one gall solution is possible for a non-trivial cc
		if ( scalar(@{$gallableTrees_ref}) == 1 ) {
			my $onlyGallableSolution = $gallableTrees_ref->[0];
			my $deletedSeq_ref = $onlyGallableSolution->[1];

			my $ccGallNode_ref = $gallNodesHash_ref->{$ccNum};
			directEdge($ccGallNode_ref, $ccNum, $deletedSeq_ref);
		}
	}
}


sub getGallStructures {
	return unless (@_ == 2);

	my $ccMatrix_ref = $_[0];
	my $cc_ref = $_[1];

	# get unique rows matrix
	my $uniqueRowsCCMatrix_ref = multiMethods->getUniqueRows($ccMatrix_ref);
	$ccMatrix_ref = multiMethods->getUniqueRowsMatrix($ccMatrix_ref, $uniqueRowsCCMatrix_ref);

 	#print "CC restricted unique rows matrix is\n";
	#for (my $i = 0; $i < scalar(@{$ccMatrix_ref}); $i++) {
	#	print "@{$ccMatrix_ref->[$i]}\n";
	#}

	my @gallableTrees;

	my $num_rows = @{$ccMatrix_ref};
	for (my $i = 0; $i < $num_rows; $i++) {
		my $copyMatrix_ref = multiMethods->getShallowMatrixCopy($ccMatrix_ref); # returns a shallow-copied matrix
		my $deletedSeq_ref = $copyMatrix_ref->[$i]; # remember sequence to be deleted

		splice @{$copyMatrix_ref}, $i, 1; # delete the i_th row from the CC restricted unique rows matrix

		#print "After deleting row ", $i + 1, "\n";
		#for (my $i = 0; $i < scalar(@{$copyMatrix_ref}); $i++) {
		#	print "@{$copyMatrix_ref->[$i]}\n";
		#}

		#  relabel first leaf as root for this algorithm to work.
		my $relabelMatrix_ref = relabelMatrix($copyMatrix_ref, 0); # relabels a deed-copied matrix using first leaf

		#print "After relabelling ...\n";
		#for (my $i = 0; $i < scalar(@{$relabelMatrix_ref}); $i++) {
		#	print "@{$relabelMatrix_ref->[$i]}\n";
		#}
		#print "\n";

		my $root = multiMethods->buildPerfPhylogeny($relabelMatrix_ref);

		if ($root == -1) { # no perfect phylogeny after removing that node
			#print "No luck with row ", $i + 1, "\n";
			next;
		}

		#print "Deleted sequence is @{$deletedSeq_ref}\n\n";
		#$root->showSubtree();
		#print "CC = @{$cc_ref}\n";

		$root->reLabelMultipleSites($cc_ref);
		$root->setZeroOneStates($copyMatrix_ref);
		#$root->showSubtree();
		#print "\n";

		my $crossOverPts_ref;
		($root, $crossOverPts_ref) = isGallableTree($root, $deletedSeq_ref, $cc_ref);

		if ( $root != 0 ) {  # tree is gallable (possibly modified tree)
			# remember the gallable tree, the corresponding deleted sequence (recombination sequence) and the crossover points
			push @gallableTrees, [$root, $deletedSeq_ref, $crossOverPts_ref];
			#$root->showSubtree();
			#print "recomb seq. = @{$deletedSeq_ref}\n";
			#print "crossovers = @{$crossOverPts_ref}\n";
		}
	}

	# remember only the two solutions with the least number of crossovers (if there are at least two solutions)
	@gallableTrees = sort { scalar(@{$b->[2]}) <=> scalar(@{$a->[2]}) } @gallableTrees;
	@gallableTrees = @gallableTrees[0..1] if ( scalar(@gallableTrees) > 1 );

	return \@gallableTrees;
}

# checks if the undirected tree is a single list and if the end nodes can be combined to get deleted sequence
# if yes, returns the $root of that (possibly modified) tree. Else, returns 0
sub isGallableTree {
	return unless (@_ == 3);

	my $root; my $deletedSeq_ref; my $cc_ref;
	($root, $deletedSeq_ref, $cc_ref) = @_;

	# check if no branching in tree except root (which can have two paths).

	my $rootChildren_ref = $root->children();
	my $numRootChildren = scalar(@{$rootChildren_ref});

	return (0,0) if ( $numRootChildren > 2 );
	croak "Tree with 3+ leaves must have children below root. Implementation bug." if ( $numRootChildren == 0 );

	my $leftEndNode_ref; my $rightEndNode_ref; my $children_ref;

	if ( $numRootChildren == 1) {
		$leftEndNode_ref = $root;
		$rightEndNode_ref = $rootChildren_ref->[0];
	}
	else { # root has exactly two branching paths; $numRootChildren = 2 here
		$leftEndNode_ref = $rootChildren_ref->[0];
		$rightEndNode_ref = $rootChildren_ref->[1];

		while (  scalar(@{$leftEndNode_ref->children()}) == 1 ) {
			$children_ref = $leftEndNode_ref->children();
			$leftEndNode_ref = $children_ref->[0];
		}

		return (0,0) if ( scalar(@{$leftEndNode_ref->children()}) > 1 ); # found a branching node other than root
	}

	while (  scalar(@{$rightEndNode_ref->children()}) == 1 ) {
		$children_ref = $rightEndNode_ref->children();
		$rightEndNode_ref = $children_ref->[0];
	}

	return (0,0) if ( scalar(@{$rightEndNode_ref->children()}) > 1 ); # found a branching node other than root

	croak "Left end leaf in tree must have sequence\n" if (not defined $leftEndNode_ref->sequence());
	croak "Right end leaf in tree must have sequence\n" if (not defined $rightEndNode_ref->sequence());

	#my $leftLeafID = $leftEndNode_ref->sequence();
	#my $rightLeafID = $rightEndNode_ref->sequence();
	#print "Left node is ", $leftLeafID + 1, "\n";
	#print "Right node is ", $rightLeafID + 1, "\n\n";

	my $leftLeafSeq_ref = $leftEndNode_ref->zeroOneState();
	my $rightLeafSeq_ref = $rightEndNode_ref->zeroOneState();

	#print "Left sequence is ", @{$leftLeafSeq_ref}, "\n";
	#print "Right sequence is ", @{$rightLeafSeq_ref}, "\n";
	#print "Deleted sequence is ", @{$deletedSeq_ref}, "\n";
	#print "Recomb point is ", $rpoint + 1, "\n\n";

	my $isGallable = 0;

	my $isLeftRight = 1;  # left end node is prefix, right end node is suffix
	$isLeftRight = 0 if ( $leftLeafSeq_ref->[0] != $deletedSeq_ref->[0] );

	my @crossOverPts;

	if ($isLeftRight == 1) {
		# assumes left as prefix and right as suffix
		getCrossOverPts($leftLeafSeq_ref, $rightLeafSeq_ref, $deletedSeq_ref, \@crossOverPts);
		$isLeftRight = 0 if ( scalar(@crossOverPts) == 0 );
	}


	if ( $isLeftRight == 1 ) {
		#print "Left prefix, Right suffix. Galled solution.\n";
		$isGallable = 1;
	}
	else {
		my $isRightLeft = 1;  # right end node is prefix, right end node is suffix
		$isRightLeft = 0 if ( $rightLeafSeq_ref->[0] != $deletedSeq_ref->[0] );

		if ($isRightLeft == 1) {
			# assumes right as prefix and left as suffix
			getCrossOverPts($rightLeafSeq_ref, $leftLeafSeq_ref, $deletedSeq_ref, \@crossOverPts);
			$isRightLeft = 0 if ( scalar(@crossOverPts) == 0 );
		}

		if ( $isRightLeft == 1 ) {
			#print "Right prefix, Left suffix. Galled solution.\n";
			$isGallable = 1;
		}
	}

	return (0,0) if ($isGallable == 0);

	if ( ($numRootChildren == 2) || ($isLeftRight == 0) ) { # redirect edges in tree so that all nodes have only one child.
				     				# For later convenience also make the prefix node made root
		my $first;
		if ($isLeftRight == 1) { # left end node is prefix, right end node is suffix
			$first = $leftEndNode_ref;
		}
		else {
			$first = $rightEndNode_ref;
		}

		my $firstAllSites = $first->allsites();
		my $firstSite = $first->site();
		$first->site(-1); # making the end leaf as root
		$first->allsites([]);
		$first->children([$first->parent()]);
		$first->parent(undef);
		my $next = $first->children()->[0];

		while ( defined $next->parent() ) {
			my $nextAllSites = $next->allsites();
			my $nextSite = $next->site();
			$next->allsites($firstAllSites);
			$next->site($firstSite);
			$next->children([$next->parent()]);
			$next->parent($first);
			$first = $next;
			$next = $first->children()->[0];
			$firstAllSites = $nextAllSites;
			$firstSite = $nextSite;
		}

		$next->allsites($firstAllSites);
		$next->site($firstSite);
		$next->parent($first);

		if ($isLeftRight == 1) { # left end prefix, right end suffix
			if ($numRootChildren == 2) {
				$next->children([$next->children()->[1]]); # only children to the right of root continue as children
			}
			else {
				$next->children([]); # no children
			}

			$root = $leftEndNode_ref;
		}
		else  { # right end prefix, left end suffix
			if ($numRootChildren == 2) {
				$next->children([$next->children()->[0]]);  # only children to the left of root continue as children
			}
			else {
				$next->children([]);
			}

			$root = $rightEndNode_ref;
		}

  		#print "Reshaped tree...\n";
		#$root->showSubtree();
		#print "\n";
	}

	# remap the crossover points in the cc-resticted matrix to those in the original matrix
	my @remapCrossOverPts;
	foreach my $point (@crossOverPts) {
		push @remapCrossOverPts, $cc_ref->[$point];
	}

	return ($root, \@remapCrossOverPts);
}


# populates the cross-over-points array if a solution exists, else makes the array empty
sub getCrossOverPts {
	return unless (@_ == 4);

	my $leftLeafSeq_ref; my $rightLeafSeq_ref; my $deletedSeq_ref; my $crossOverPts_ref;
	($leftLeafSeq_ref, $rightLeafSeq_ref, $deletedSeq_ref, $crossOverPts_ref) = @_;

	my $mode = 0; # look at left (prefix) initally

	for (my $i = 1; $i < scalar(@{$deletedSeq_ref}); $i++) {

		if ($mode == 0) {
			if ( $leftLeafSeq_ref->[$i] != $deletedSeq_ref->[$i] ) {
				my $lastCrossOverID = scalar(@{$crossOverPts_ref}) - 1;
				if ( ($lastCrossOverID >= 0) && ($crossOverPts_ref->[$lastCrossOverID] == $i) ) { # both left and right sequences mismatch
					@{$crossOverPts_ref} = ();
					#print "1. Left not prefix. Mismatch found.\n";
					last;
				}

				push @{$crossOverPts_ref}, $i;
				$i--;
				$mode = 1; # look at right (suffix) the next time
			}
		}
		else {
			if ( $rightLeafSeq_ref->[$i] != $deletedSeq_ref->[$i] ) {
				my $lastCrossOverID = scalar(@{$crossOverPts_ref}) - 1;
				if ( ($lastCrossOverID >= 0) && ($crossOverPts_ref->[$lastCrossOverID] == $i) ) { # both left and right sequences mismatch
					@{$crossOverPts_ref} = ();
					#print "2. Left not prefix. Mismatch found.\n";
					last;
				}

				push @{$crossOverPts_ref}, $i;
				$i--;
				$mode = 0;
			}
		}
	}

}


# relabel matrix with $leafid as root, ie, columns with a 1 in that sequence get changed so that 0s->1s and 1s->0s
sub relabelMatrix {
	return unless (@_ == 2);

	my $given_ref = $_[0];
	my $leaf_id = $_[1];

	my $return_ref = multiMethods->getDeepMatrixCopy($given_ref);
	my $rowCount = @{$given_ref};
	my $colCount = @{$given_ref->[0]};

	for (my $j = 0; $j < $colCount; $j++) {

		my $site = ${$given_ref->[$leaf_id]}[$j];

		if ($site == 1) { # relabel the column j of matrix
			for (my $i = 0; $i < $rowCount; $i++) {
				${$return_ref->[$i]}[$j] = 1 - ${$return_ref->[$i]}[$j]; # 0->1 and 1->0
			}
		}
	}

	return $return_ref;
}


# get a reference to an array of CC restricted matrices
sub getCCMatrices
{
	return unless @_ == 2;
	my $matrix_ref; my $ccs_ref;

	($matrix_ref, $ccs_ref) = @_;

	my $num_cc = @{$ccs_ref};
	my $ccMatrices_ref;

	#print "Computing matrices restricted to the columns in each CC.\n";
	for (my $i = 0; $i < $num_cc; $i++) {
		push @{$ccMatrices_ref}, multiMethods->getCCMatrix($matrix_ref, \@{$ccs_ref->[$i]});
	}

	return $ccMatrices_ref;
}


sub getSplitMatrix
{
	return unless @_ == 1;
	my $ccMatrices_ref = $_[0];

	my $num_cc = @{$ccMatrices_ref};
	my $num_rows = @{$ccMatrices_ref->[0]};
	my @splitMatrix; my @splitLeaves; my @splitMatrixColumnMap;

	for (my $i = 0; $i < $num_cc; $i++) {

		#print "The restricted matrix for CC ", ($i + 1) ," is\n";
		#for (my $j = 0; $j < $num_rows; $j++) {
		#	print @{$ccMatrices_ref->[$i][$j]}, " (row ", $j+1, ")", "\n";
		#}
		#print "\n";

		my $uniqueCCRows_ref = multiMethods->groupIdenticalRows($ccMatrices_ref->[$i]);

		# add columns in split matrix corresponding to each CC
		my @sortedKeys = sort { $a <=> $b } (keys %{$uniqueCCRows_ref} );

		if ( scalar(@sortedKeys) > 2 ) { # non-trivial CC
			foreach my $item ( @sortedKeys ) {
				#print "$item --maps to rows--> @{$uniqueCCRows_ref->{$item}} \n";

				# remember those splits that split exactly one leaf
				my $numSplitLeaves = @{$uniqueCCRows_ref->{$item}};
				if ($numSplitLeaves == 1) {
					push @splitLeaves, ${$uniqueCCRows_ref->{$item}}[0];
				}

				my $rowNumCCMatrix = ${$uniqueCCRows_ref->{$item}}[0];
				#print "CC # = $i, seq =  @{$ccMatrices_ref->[$i][$rowNumCCMatrix]}\n";
				push @splitMatrixColumnMap, [$i, \@{$ccMatrices_ref->[$i][$rowNumCCMatrix]}];
				multiMethods->addColumn(\@splitMatrix, \@{$uniqueCCRows_ref->{$item}}, $num_rows);
			}
		}
		elsif ( scalar(@sortedKeys) == 2 ) { # trivial CC defining a split twice. So, adding only one column to split matrix

			#print "0 --> @{$uniqueCCRows_ref->{0}} \n";
			#print "1 --> @{$uniqueCCRows_ref->{1}} \n";

			# remember those splits that split exactly one leaf
			my $numSplitLeaves = @{$uniqueCCRows_ref->{$sortedKeys[0]}};

			# either the unique row 0 or the unique row 1 can split exactly one leaf
			if ($numSplitLeaves == 1) {
				push @splitLeaves, ${$uniqueCCRows_ref->{$sortedKeys[0]}}[0];
			}
			elsif ($numSplitLeaves == $num_rows - 1) {
				push @splitLeaves, ${$uniqueCCRows_ref->{$sortedKeys[1]}}[0];
			}

			my $rowNumCCMatrix = ${$uniqueCCRows_ref->{$sortedKeys[0]}}[0];
			#print "CC # = $i, seq = @{$ccMatrices_ref->[$i][$rowNumCCMatrix]}\n";
			push @splitMatrixColumnMap, [$i, \@{$ccMatrices_ref->[$i][$rowNumCCMatrix]}];
			multiMethods->addColumn(\@splitMatrix, \@{$uniqueCCRows_ref->{$sortedKeys[0]}}, $num_rows);
		}
		else {
			croak "A column with all 0s or all 1s. Must not have occured\n";
		}

		#print "\n\n";
	}
	#print "\n";


	#print "The matrix of splits (possibly without all leaf splits) is\n";
	#for (my $j = 0; $j < $num_rows; $j++) {
	#	print @{$splitMatrix[$j]}, " (row ", $j+1, ")", "\n";
	#}
	#print "\n";

	#print "The split leaves are @splitLeaves\n";
	#print "Leaves to be split are\n";

	# add a column in the split matrix for each leaf that is not yet split
	@splitLeaves = (sort { $a <=> $b } @splitLeaves);
	my $rowID = 0; my $numSplitLeaves = @splitLeaves;
	for (my $i = 0; $i < $num_rows; $i++) {
		if ( ($rowID < $numSplitLeaves)  && ($splitLeaves[$rowID] != $i) ) {
			# leaf $i must be split
			#print "$i ";
			#print "CC # = -1, leaf # = $i\n";
			push @splitMatrixColumnMap, [-1, [$i]];
			multiMethods->addColumn(\@splitMatrix, [$i], $num_rows);
		}
		else {
			$rowID++;
			while ( ($splitLeaves[$rowID] == $splitLeaves[$rowID-1]) && ($rowID < $numSplitLeaves) ) { # go to next unique leaf
				$rowID++;
			}
		}
	}
	#print "\n\n";

	return (\@splitMatrix, \@splitMatrixColumnMap);
}


# returns 1 if the array elements are equal, 0 otherwise
sub isEqual
{
	my $array1_ref = $_[0];
	my $array2_ref = $_[1];

	my $num1_rows = @{$array1_ref};
	my $num2_rows = @{$array2_ref};

	if ($num1_rows != $num2_rows) {
		#print "Arrays unequal - # of elements different\n";
		return 0;
	}

	for (my $i = 0; $i < $num1_rows; $i++) {
		if ( $array1_ref->[$i] != $array2_ref->[$i] ) {
			#print "Arrays unequal - differing at element $i\n";
			return 0;
		}
	}

	# arrays are equal
	return 1;
}


# Returns the reference to an array of row numbers(species) with both the (characteristics)column values as 1
sub getIntersect
{
	my $matrix_ref = $_[0];
	my $col1 = $_[1];
	my $col2 = $_[2];

	print "getIntersect - columns $col1 and $col2\n";

	my $num_rows = @{$matrix_ref};
	my @species12;

        for (my $i = 0; $i < $num_rows; $i++) {
		#print "'$matrix_ref->[$i][$col1], $matrix_ref->[$i][$col2]'\n";

		if ( ($matrix_ref->[$i][$col1] == 1) &&
		     ($matrix_ref->[$i][$col2] == 1) ) {
			push @species12, $i;
		}
	}

	return \@species12;
}


# Read file into a matrix of matrices
sub matrix_read_file
{
	my ($filename) = $_[0];
	open (FILE, $filename) || die "Could not open $filename";
	while ( defined( $line = <FILE> ) )
	{
		chomp($line);
		next if $line =~ /^\s*$/; #skip blank lines
		if ($line =~ /^([A-Za-z]\w*)/) {
			# assuming matrix names begin with letters
			$matrix_name = $line;
			push @matrices, $matrix_name; # insert all matrix names
		}
		else {
			my $texlen = length($line);
			my $newtext = "";
			for (my $i = 0; $i < $texlen; $i++) {
				$char = substr($line, $i,1);
				next if ($char !~ /[0-1]/);
				$newtext .= $char;
				$newtext .= " ";
			}

			my (@row) = split (/\s+/, $newtext);
			push @{$matrix_name}, \@row; #insert row arrays
		}
	}
	close(FILE);

	return (\@matrices);
}


# Test for case (a)
# Load the conflicting columns for each column into a hash of matrices
sub findConflicts
{
	my $matrix_ref = $_[0];
        #print "Testing all column pairs for conflict\n";

	my $num_rows = @{$matrix_ref};
	my $num_cols = @{$matrix_ref->[0]}; # assuming all rows have same
					    # number of columns

	# Perfect phylogenetic network without recombination
	# exists only if none of the columns conflict
	my $bNoConflicts = 1;

	# Initializing a hash of anonymous empty arrays
	my %conflictColumns;
	for (my $i = 0; $i < $num_cols; $i++) {
		$conflictColumns{$i} = [];
	}

	for (my $i = 0; $i < $num_cols; $i++) {
		for ($j = $i + 1; $j < $num_cols; $j++) {
			if ( isConflictCols($matrix_ref, $i, $j) ) {
				$bNoConflicts = 0;

				# remember the set of columns that conflict with each column
				push @{$conflictColumns{$i}}, $j;
				push @{$conflictColumns{$j}}, $i;
			}
		}
	}

	print "\n";
	return ($bNoConflicts, \%conflictColumns);
}


# Tests if for a given pair of columns, the value pairs (0, 0) (0,1) (1,0) and (1,1) exist.
# If so, returns 1 else returns 0.
sub isConflictCols
{
	my $matrix_ref = $_[0];
	my $col1 = $_[1];
	my $col2 = $_[2];

	#print "Testing columns $col1 and $col2\n";

	my $num_rows = @{$matrix_ref};
	my $oneone; my $zeroone; my $onezero; my $zerozero;
	$oneone = $zeroone = $onezero = $zerozero = 0;

        for (my $i = 0; $i < $num_rows; $i++) {
		if ( ($matrix_ref->[$i][$col1] == 1) &&
		     ($matrix_ref->[$i][$col2] == 1) ) {
			$oneone = 1;
		}
		elsif ( ($matrix_ref->[$i][$col1] == 0) &&
                          ($matrix_ref->[$i][$col2] == 1) ) {
                        $zeroone = 1;
                }
		elsif ( ($matrix_ref->[$i][$col1] == 1) &&
                          ($matrix_ref->[$i][$col2] == 0) ) {
                        $onezero = 1;
                }
		else {
			$zerozero = 1;
		}

		if ($oneone && $zeroone && $onezero && $zerozero) {
			return 1;
		}
	}

	return 0;
}

sub getConnectedComponents
{
	my $matrix_ref = $_[0];
	my $confColHash_ref = $_[1];

	my @ccHelper; #used to compute the connected components

	my $num_cols =  @{$matrix_ref->[0]};
	for (my $i = 0; $i < $num_cols; $i++) {
		push @ccHelper, $i;
	}

	#print "Initial state of connected component helper\n";
	#print "@ccHelper\n";

	for (my $col1 = 0; $col1 < $num_cols; $col1++) {
		my $num_col1 = @{$confColHash_ref->{$col1}}; # number of columns conflicting with column 'col1'

		for (my $i = 0; $i < $num_col1; $i++) {
			my $col2 = ${$confColHash_ref->{$col1}}[$i];

			# (col1, col2) is an edge in the conflict graph
			if ($ccHelper[$col1] != $ccHelper[$col2])
			{
				my $min; my $max;
				if ($ccHelper[$col1] < $ccHelper[$col2])
				{
					$min = $ccHelper[$col1];
					$max = $ccHelper[$col2]
				}
				else
				{
					$max = $ccHelper[$col1];
					$min = $ccHelper[$col2]
				}

				# replace all occurences of $max with $min
				for (my $j = 0; $j < $num_cols; $j++) {
					if ($ccHelper[$j] == $max) {
						$ccHelper[$j] = $min;
					}
				}
			}
		}
	}

	#print "Final state of connected component helper\n";
	#print "@ccHelper\n";

	#create an array of connected components
	my @ccArray;

	for (my $i = 0; $i < $num_cols; $i++) {
		if ($ccHelper[$i] == $i) {
			my @row;
			for (my $j = 0; $j < $num_cols; $j++) {
				if ($ccHelper[$j] == $i) {
					push @row, $j;
				}
			}
			push @ccArray, \@row;
		}
	}

	#print "The connected components are below\n";
	#my $num_cc = @ccArray;
	#for (my $i = 0; $i < $num_cc; $i++) {
	#	print "@{$ccArray[$i]}\n";
	#}

	return \@ccArray;
}


# parameters(original matrix, unique columns array)
#returns a matrix that retains only those unique columns of the original matrix
sub getUniqueMatrix
{
	my $original_matrix_ref = $_[0];
	my %uniqueCols = %{$_[1]};
	my %uniqueRows = %{$_[2]};

	my @uniqueMatrix;

	foreach my $rowItem (sort { $a <=> $b } (keys % uniqueRows) ) {
		my $newtext = "";

		foreach my $colItem (sort { $a <=> $b } (keys % uniqueCols) ) {
			$newtext .= $original_matrix_ref->[$rowItem][$colItem];
			$newtext .= " ";
		}

		my (@row) = split (/\s+/, $newtext);
		push @uniqueMatrix, \@row;
	}

	#my $num_rows = @uniqueMatrix;
	# printing each row of the matrix
        #for (my $i = 0; $i < $num_rows; $i++) {
	#	print "@{$uniqueMatrix[$i]}\n";
	#}
	#print "\n";

	return \@uniqueMatrix;
}
