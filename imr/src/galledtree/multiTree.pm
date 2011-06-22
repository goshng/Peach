package multiTree::Node;
use Carp;


sub new {
	my $type = shift;
	my $class = ref($type) || $type;
	croak "Node::new needs at least 1 argument." unless @_ > 0;
	my $site = shift;

	my $self = {};
	$self->{site} = $site;
	$self->{allsites} = []; # sites corresponding to node. site + duplicates = allsites
	$self->{sequence} = undef;
	$self->{parent} = undef;
	$self->{seekParent} = 0; # states whether or not to direct the edge to parent. 0 is undirected, 1 is directed.
	$self->{children} = [];
	$self->{seekChildren} = {};  # states whether or not to direct the edges to each chid. 0 is undirected, 1 is directed.
	$self->{onsites} = []; # sites that will be 1 at this node (all parent sites propagated will also be 1 here)
	$self->{zeroOneState} = []; # stores 011000xxx.. the state of 0s and 1s in node

	bless $self, $class;

	return $self;
}

sub seekParent {
	my $self = shift;
	$self->{seekParent} = shift if @_;
	return $self->{seekParent};
}

sub seekChildren {
	my $self = shift;
	return $self->{seekChildren};
}

# sequence of recombination node
sub zeroOneState {
	my $self = shift;
	$self->{zeroOneState} = shift if @_;

	return $self->{zeroOneState};
}

sub site {
	my $self = shift;
	$self->{site} = shift if @_;
	return $self->{site};
}

sub sequence {
	my $self = shift;
	$self->{sequence} = shift if @_;
	return $self->{sequence};
}

sub parent {
	my $self = shift;
	$self->{parent} = shift if @_;
	return $self->{parent};
}

sub children {
	my $self = shift;
	$self->{children} = shift if @_;
	return $self->{children};
}

sub onsites {
	my $self = shift;
	$self->{onsites} = shift if @_;
	return $self->{onsites};
}

sub allsites {
	my $self = shift;
	$self->{allsites} = shift if @_;
	return $self->{allsites};
}

# get the list of 'on' sites of parent. Add my 'on' site to that list.
sub setOnsites {
	my $self = shift;
	return unless @_ == 4;

	my $num_cols; my $parentOns_ref; my $matrix_ref; my $rootSeqID;
	($num_cols, $parentOns_ref, $matrix_ref, $rootSeqID) = @_;

	if ( $self->site() >= 0 ) {
		my $allsites_ref = $self->allsites();
		my @onsites = @{[@{$parentOns_ref}, @{$allsites_ref}]};

		$self->onsites(\@onsites); # sites that are 'on' at this node
	}
	else { # avoid propagating site -1 of 'root' node
		$self->onsites($parentOns_ref);
	}

	#print "On sites of node ", $self->site(), " are @{$self->onsites()}\n";

	if ( defined $self->sequence() ) {
		#print "Checking equality for site : ", $self->site() + 1, ", on sites = @{$self->{onsites}}\n";
		# check if the set of 'on' sites matches with the corresponding input sequence
		if ( multiMethods->checkSequenceEquality($matrix_ref->[$rootSeqID], $self->onsites(), $matrix_ref->[$self->sequence()]) == 0 ) {
			print "Failed corresponding sequence equality for site ", $self->site() + 1, "\n";
			return -1; # if not equal
		}

		$self->zeroOneState($matrix_ref->[$self->sequence()]);
	}
	else {
		my $zeroOneState_ref = multiMethods->getZeroOneState($matrix_ref->[$rootSeqID], $num_cols, $self->onsites());
		$self->zeroOneState($zeroOneState_ref);
	}
}


sub findOnSites {
	my $self = shift;
	return unless @_ == 6;

	my $num_cols; my $parentOns_ref; my $matrix_ref; my $rootSeqID; my $site_counter_ref; my $leaf_counter_ref;
	($num_cols, $parentOns_ref, $matrix_ref, $rootSeqID, $site_counter_ref, $leaf_counter_ref) = @_;

	if ( $self->site() >= 0 ) { # increment the counter for these sites
		my @allsites = @{$self->allsites()};
		for (my $i = 0; $i < scalar(@allsites); $i++) {
			$site_counter_ref->[$allsites[$i]]++;
		}
	}

	if ( defined $self->sequence() ) { # increment the counter for this leaf sequence
		$leaf_counter_ref->[$self->sequence()]++;
	}

	if ( $self->setOnsites($num_cols, $parentOns_ref, $matrix_ref, $rootSeqID) == -1) {
		return -1;
	}

	foreach my $child ( @{$self->children()} ) {
		if ( $child->findOnSites($num_cols, \@{$self->onsites()}, $matrix_ref,
						$rootSeqID, $site_counter_ref, $leaf_counter_ref) == -1 ) {
			return -1;
		}
	}

}


# adds only one child, even if multiple children are presented
sub addChild {
	my $self = shift;
	if ( @_ > 0 ) {
		my $child = shift;

		#print "Adding ", $child->site(), " as child of ", $self->site(), "\n";
		#print "Adding ...\n";
		#$child->display();
		#print " as child of ...\n";
		#$self->display(); print " - $self\n",

		# Newly created edges are bidirectional
		push @{$self->children()}, $child;
		my $hash_ref = $self->seekChildren();
		$hash_ref->{$child} = 1;  # allow parent->child access
		$child->parent($self);
		$child->seekParent(1); # allow child->parent access
	}
}


# some subtrees of the perfect phylogeny displaying site mutations have no rows
# associated with them. Remove those extraneous subtrees.
# NOTE: UNNECESSARILY HEAVY-WEIGHT IMPLEMENTATION. SITES WITH NO SEQUENCES CANNOT
# FORM A NETWORK. SO, THEY WILL ALL BE CHILDREN OF ROOT. ENOUGH TO EXAMINE JUST
# THEM AND REMOVE.
sub removeSequencelessSubTree {
	my $self = shift;
	my $childIndexOfParent = shift;

	my $size = @{$self->children()};
	my @nodes = @{$self->children()};
	for (my $i = 0; $i < @{$self->children()}; $i++) {
		my $node = $nodes[$i];
		if ( $node->removeSequencelessSubTree($i) == 1 ) {
			$i--;
		}
	}

	my $deleted = 0;

	if ( (defined $self->sequence()) || (@{$self->children()} > 0) ) {
		# do not delete if sequence found or child exists
		#print "Not deleting ..."; $self->display();
	}
	elsif ( $childIndexOfParent >= 0 ) { # Do not remove the only root with no parent.
		# disconnect from parent.

		my $parent = $self->parent(); #print "Deleting child # - $childIndexOfParent of "; $parent->display();
		my $children = $parent->children();
		splice @{$children}, $childIndexOfParent, 1;

		$deleted = 1;
	}

	return $deleted;
}

# increasing sites and sequence numbers by 1 to be 1-based instead of being traditionally 0-based
sub showSubtree {
	my $self = shift;
	my $indent = shift || 0;
	my $I = ' ' x $indent;
	print $I;

	if ($self->site() == -1) {  # check if site is root
		print "root";
	}
	else {
		my @allsites = @{$self->allsites()};
		if ( scalar(@allsites) > 0 ) {
			print $allsites[0] + 1;

			for (my $i = 1; $i < scalar(@allsites); $i++) {
	 			print ",", $allsites[$i] + 1;
	 		}
		}
		else {
			print "Incorrect: no sites are 'on'\n";
		}
	 }

	print " (row ";
	(defined $self->sequence()) ? print $self->sequence() + 1 : print "?"; # check if sequence exists
	print ") ", @{$self->zeroOneState()}, "\n";

	my $size = @{$self->children()};
	my @nodes = @{$self->children()};
	for (my $i = 0; $i < $size; $i++) {
		my $node = $nodes[$i];
		$node->showSubtree($indent+2);
	}
}


# every edge label is of the form ( a, b[] ).
sub showNodeCompoundEdgeLabels {
	my $self = shift;

	if ($self->site() == -1) {  # check if site is root
		print "root";
	}
	else {
		my @allsites = @{$self->allsites()};
		if ( scalar(@allsites) > 0 ) {
			my $item = $allsites[0];
			if (defined $item->[1]) {
				print "($item->[0], @{$item->[1]})";
			}
			else {
				print $allsites[0] + 1;
			}

			for (my $i = 1; $i < scalar(@allsites); $i++) {
				$item = $allsites[$i];
				if (defined $item->[1]) {
					print ",($item->[0], @{$item->[1]})";
				}
				else {
					print ",", $allsites[$i] + 1;
				}
	 		}
		}
		else {
			print "Incorrect: no sites are 'on'\n";
		}
	 }

	print " (row ";
	(defined $self->sequence()) ? print $self->sequence() + 1 : print "?"; # check if sequence exists
	print ") ", @{$self->zeroOneState()}, "\n";
}


# returns 1 if the phylogenetic network has a gall with multiple crossovers; 0 otherwise
sub hasMultipleCrossovers {
	my $self = shift;

	my @nodes = @{$self->children()};
	foreach my $node ( @nodes ) {
		return 1 if ($node->hasMultipleCrossovers() == 1);
	}

	return 0;
}



# every edge label is of the form ( a, b[] ).
sub showSubtreeCompoundEdgeLabels {
	my $self = shift;
	my $indent = shift || 0;
	my $I = ' ' x $indent;
	print $I;

	if ($self->site() == -1) {  # check if site is root
		print "root";
	}
	else {
		my @allsites = @{$self->allsites()};
		if ( scalar(@allsites) > 0 ) {
			my $item = $allsites[0];
			if (defined $item->[1]) {
				print "($item->[0], @{$item->[1]})";
			}
			else {
				print $allsites[0] + 1;
			}

			for (my $i = 1; $i < scalar(@allsites); $i++) {
				$item = $allsites[$i];
				if (defined $item->[1]) {
					print ",($item->[0], @{$item->[1]})";
				}
				else {
					print ",", $allsites[$i] + 1;
				}
	 		}
		}
		else {
			print "Incorrect: no sites are 'on'\n";
		}
	 }

	print " (row ";
	(defined $self->sequence()) ? print $self->sequence() + 1 : print "?"; # check if sequence exists
	print ") ", @{$self->zeroOneState()}, "\n";

	my $size = @{$self->children()};
	my @nodes = @{$self->children()};
	for (my $i = 0; $i < $size; $i++) {
		my $node = $nodes[$i];
		$node->showSubtreeCompoundEdgeLabels($indent+2);
	}
}


# increasing sites and sequence numbers by 1 to be 1-based instead of being traditionally 0-based
sub displaySubtreeSingleEdgeLabel {
	my $self = shift;
	my $indent = shift || 0;
	my $I = ' ' x $indent;
	print $I;

	if ($self->site() == -1) {  # check if site is root
		print "root";
	}
	else {
		print $self->site() + 1;
	}

	print " (row ";
	(defined $self->sequence()) ? print $self->sequence() + 1 : print "?"; # check if sequence exists
	print ") ", @{$self->zeroOneState()}, "\n";

	my $size = @{$self->children()};
	my @nodes = @{$self->children()};
	for (my $i = 0; $i < $size; $i++) {
		my $node = $nodes[$i];
		$node->displaySubtree($indent+2);
	}
}


# add leaves corresponding to nodes in the subtree(with the node as root)
sub addSubtreeLeaves
{
	my $self = shift;
	my $hash_ref = shift if @_;

	return unless defined ($hash_ref);

	$site = $self->site();
	if ( exists $hash_ref->{$site} ) {
		$self->sequence($hash_ref->{$site});
	}

	my $size = @{$self->children()};
	my @childNodes = @{$self->children()};
	for (my $i = 0; $i < $size; $i++) {
		my $child = $childNodes[$i];
		$child->addSubtreeLeaves($hash_ref);
	}
}


# Algorithm built perfect phylogeny for column-reordered matrix M'.
# So, relabelling the nodes to correspond to the original matrix M.
sub reLabelSites {
	my $self = shift;
	my $order_ref = shift if @_;

	return unless defined ($order_ref);

	my $oldSite = $self->site();
	if ( ($oldSite >= 0) && ($oldSite < scalar(@{$order_ref})) ) {
		$self->site($order_ref->[$oldSite]);
	}

	my $size = @{$self->children()};
	my @childNodes = @{$self->children()};
	for (my $i = 0; $i < $size; $i++) {
		my $child = $childNodes[$i];
		$child->reLabelSites($order_ref);
	}
}


# relabels multiple labels on edges in T_bar to CC based labels
sub reLabelMultipleSites {
	my $self = shift;
	my $order_ref = shift if @_;

	return unless defined ($order_ref);

	my $oldSite = $self->site();
	if ($oldSite != -1) {
		$self->site($order_ref->[$oldSite]);
	}

	my @allsites = @{$self->allsites()};
	if ( scalar(@allsites) > 0 ) {
		for (my $i = 0; $i < scalar(@allsites); $i++) {
			my $oldSite = $allsites[$i];
			if ( ($oldSite >= 0) && ($oldSite < scalar(@{$order_ref})) ) {
				#print "$allsites[$i]-->($order_ref->[$oldSite][0], @{$order_ref->[$oldSite][1]})\n";
				$allsites[$i] = $order_ref->[$oldSite];
			}
		}
	}
	$self->allsites(\@allsites);

	my $size = @{$self->children()};
	my @childNodes = @{$self->children()};
	for (my $i = 0; $i < $size; $i++) {
		my $child = $childNodes[$i];
		$child->reLabelMultipleSites($order_ref);
	}
}


sub reLabelSequences {
	my $self = shift;
	my $order_ref = shift if @_;

	return unless defined ($order_ref);

	my $oldSequence = $self->sequence();
	if ( defined $oldSequence ) {
		$self->sequence($order_ref->[$oldSequence]);
	}

	my $size = @{$self->children()};
	my @childNodes = @{$self->children()};
	for (my $i = 0; $i < $size; $i++) {
		my $child = $childNodes[$i];
		$child->reLabelSequences($order_ref);
	}
}


sub setZeroOneStates {
	my $self = shift;
	my $matrix_ref = shift if @_;

	return unless defined ($matrix_ref);

	my $sequence = $self->sequence();
	if ( defined $sequence ) {
		if ( ($sequence >= 0) && ($sequence < scalar(@{$matrix_ref})) ) {
			$self->zeroOneState($matrix_ref->[$sequence]);
		}
	}

	my $size = @{$self->children()};
	my @childNodes = @{$self->children()};
	for (my $i = 0; $i < $size; $i++) {
		my $child = $childNodes[$i];
		$child->setZeroOneStates($matrix_ref);
	}
}


sub addDuplicateEdgeLabels {
	my $self = shift;
	return unless @_ == 1;

	my $duplicateColsHash_ref = $_[0];

	my $site = $self->site(); # get new site
	my @array;

	push @array, $site; # original site

	if ( defined $duplicateColsHash_ref->{$site} ) { # if site has duplicates
		push @array, @{$duplicateColsHash_ref->{$site}}; # add duplicates of the site
	}

	$self->allsites(\@array); # site + duplicates = all sites

	foreach my $child ( @{$self->children()} ) {
		$child->addDuplicateEdgeLabels($duplicateColsHash_ref);
	}
}


sub getGallNodes {
	my $self = shift;
	return unless @_ == 2;

	my $gallNodesHash_ref = $_[0];
	my $num_cc = $_[1];

	# if interior node
	if ( (scalar(@{$self->children()}) > 0) &&  (defined $self->parent()) ) {
		# test if expandable gall node

		# initialize a set of flags to 1. (assuming all labels are present)
		my @ccsFlag = (0..$num_cc-1);
		for (my $i = 0; $i < scalar(@ccsFlag); $i++) {
			$ccsFlag[$i] = 1;
		}

		# initially assume that no labels are on an edge
		my @edgeLabelsFlag = (0..$num_cc-1);
		for (my $i = 0; $i < scalar(@edgeLabelsFlag); $i++) {
			$edgeLabelsFlag[$i] = 0;
		}

		foreach my $label ( @{$self->allsites()} )  {
			my $cc_num = $label->[0];
			if ($cc_num >= 0 && $cc_num < $num_cc) { # cc label, not a leaf label (leaf labels have cc# = -1)
				$edgeLabelsFlag[$cc_num] = 1; # remember all labels on an edge
			}
		}

		# retain only those labels that are present in the touching edge leading from its parent
		for (my $i = 0; $i < $num_cc; $i++) {
			$ccsFlag[$i] = $ccsFlag[$i] * $edgeLabelsFlag[$i];
		}

		# further retain only those labels that are present in a touching edge leading to its children
		foreach my $child ( @{$self->children()} ) {

			for (my $i = 0; $i < scalar(@edgeLabelsFlag); $i++) {
				$edgeLabelsFlag[$i] = 0;
			}

			foreach my $label ( @{$child->allsites()} )  {
				my $cc_num = $label->[0];
				if ($cc_num >= 0 && $cc_num < $num_cc) { # cc label, not a leaf label
					$edgeLabelsFlag[$cc_num] = 1;
				}
			}

			# retain the 1 on only those cc flags that are in touching ccs
			for (my $i = 0; $i < $num_cc; $i++) {
				$ccsFlag[$i] = $ccsFlag[$i] * $edgeLabelsFlag[$i];
			}
		}

		# only one cc label must eventually be retained for a node
		# this must be the only node in T_bar which retains that label

		my $found = 0;
		for (my $i = 0; $i < scalar(@ccsFlag); $i++) {
			if ( $ccsFlag[$i] == 1 ) {
				if ( not defined $gallNodesHash_ref->{$i} ) {
					if ( $found == 1) {
						croak "Cannot expand a node into more than a gall. Something wrong!\n";
					}

					$found = 1;
					$gallNodesHash_ref->{$i} = $self;
				}
				else {
					croak "Cannot have more than one node expanding into the same gall. Something wrong!\n";
				}

			}

		}

	}

	# test each child node
	foreach my $child ( @{$self->children()} ) {
		$child->getGallNodes($gallNodesHash_ref, $num_cc);
	}

}


# collect all the nodes in the subtree for which it is the root
sub collectAllNodes {
	my $self = shift;
	return unless @_ == 1;

	my $collectedNodes_ref = $_[0];
	push @{$collectedNodes_ref}, $self;

	foreach my $child ( @{$self->children()} ) {
		$child->collectAllNodes($collectedNodes_ref);
	}
}


# returns 1 if node can become the root; 0 otherwise
sub isRootable {
	my $self = shift;

	return 0 if (not defined $self->sequence()); # choosing only a leaf node as root for convenience. Any galled-tree has a corresponding
						     # galled-tree with one of the leaves as root. So this restraint does not contradict the goal.
	return 0 if ( $self->canAccessAllAbove() == 0 ); # fail if cannot access all the nodes below
	return 0 if ( $self->canAccessAllBelow() == 0 ); # fail if cannot access all the nodes above
	return 1;
}

# returns 1 if node can access all nodes in subtree for which it is the root; 0 otherwise
sub canAccessAllBelow {
	my $self = shift;

	my $childSeekHash_ref = $self->seekChildren();

	my @children = @{$self->children()};
	my $numChildren = scalar(@children);

	for (my $i = 0; $i < $numChildren; $i++) {
		my $child_ref = $children[$i];

		#$child_ref->showNodeCompoundEdgeLabels();
		return 0 if ( $childSeekHash_ref->{$child_ref} == 0 ); # cannot access a child
		return 0 if ( $child_ref->canAccessAllBelow() == 0 ); # child cannot access its children
	}

	return 1; # can access all children, each child can access all its children
}

# returns 1 if the node can access all nodes except those in the subtree for which one of its children is the root; 0 otherwise
sub canAccessAllButOneBelow {
	my $self = shift;
	return unless @_ == 1;

	my $unExploredChild_ref = $_[0];
	my $childSeekHash_ref = $self->seekChildren();

	my @children = @{$self->children()};
	my $numChildren = scalar(@children);

	for (my $i = 0; $i < $numChildren; $i++) {
		my $child_ref = $children[$i];
		next if ($child_ref == $unExploredChild_ref); # don't explore path leading to this child
		#$child_ref->showNodeCompoundEdgeLabels();

		return 0 if ( $childSeekHash_ref->{$child_ref} == 0 ); # cannot access a child
		return 0 if ( $child_ref->canAccessAllBelow() == 0 ); # child cannot access its children
	}

	return 1; # can access all children, each child can access all its children
}

# returns 1 if the node can access all nodes except those that are in the subtree for which it is the root; 0 otherwise
sub canAccessAllAbove {
	my $self = shift;

	my $parent = $self->parent();
	if (defined $parent) {
		#$parent->showNodeCompoundEdgeLabels();
		return 0 if ( $self->seekParent() == 0 ); # parent not accessible
		return 0 if ( $parent->canAccessAllAbove() == 0 ); # parent cannot access all above
		return 0 if ( $parent->canAccessAllButOneBelow($self) == 0 ); # parent cannot access all its other children
	}

	return 1;
}

# directs all edges away from this node to make it as the new root
sub directEdgesAway {
	my $self = shift;

	# reverse the direction of all edges up to root
	$self->pointToParent();

	# making $self as root
	$self->parent(undef);
	$self->seekParent(0);
	$self->site(-1);
	$self->allsites([]);
}

# reverses the direction of all edges from $self up to the root
sub pointToParent {
	my $self = shift;

	my $parent = $self->parent();
	if (defined $parent) {
		# make the parent point to its ancestor (grandparent)
		$parent->pointToParent();

		# add previous parent as its child
		$self->addChild($parent);

		# directing edge away from $self
		$parent->site($self->site());
		$parent->allsites($self->allsites());

		# previous parent no longer has $self as its child
		my @children = @{$parent->children()};
		my $numChildren = scalar(@children);
		for (my $i = 0; $i < $numChildren; $i++) {
			my $child_ref = $children[$i];
			if ($child_ref == $self) {
				splice @{$parent->children()}, $i, 1;
				last;
			}
		}
	}
}


# T_bar nodes with an immediately branching leaf node can absorb that leaf into themselves
# unnecessary immediately branching child with no sequence (since parent is itself that sequence) can be deleted
sub mergeSelfLeaf {
	my $self = shift;

	my $children_ref = $self->children();
	my $numChildren = scalar(@{$children_ref});

	if (not defined $self->sequence()) {
		for (my $i = 0; $i < $numChildren; $i++) {
			my $child_ref = $children_ref->[$i];
			my @allsites = @{$child_ref->allsites()};
			if ( (scalar(@allsites) == 1) && ($allsites[0]->[0] == -1) ) { # child node is an immediately branching leaf
				croak "Immediately branching child must know its sequence. Implementation bug!" if (not defined $child_ref->sequence());
				$self->sequence($child_ref->sequence()); # absorbing leaf node into itself
				splice @{$children_ref}, $i, 1; # deleting the unnecessary child
				last;
			}
		}
	}
	else {
		for (my $i = 0; $i < $numChildren; $i++) {
			my $child_ref = $children_ref->[$i];
			my @allsites = @{$child_ref->allsites()};
			if ( (scalar(@allsites) == 1) && ($allsites[0]->[0] == -1) && (not defined $child_ref->sequence()) ) {
				 # child node is an unnecessary immediately branching leaf defining no sequence since the parent defines that sequence
				 # merge that child node with the parent and let the parent adopt the child's children (grand children) as its own children
				 foreach my $grandChild ( @{$child_ref->children()} ) {
				 	$self->addChild($grandChild); # adopting grand children
				 }
				splice @{$children_ref}, $i, 1; # deleting the unnecessary child
				last;
			}
		}
	}

	# check the rest of the tree
	$children_ref = $self->children();
	$numChildren = scalar(@{$children_ref});
	for (my $i = 0; $i < $numChildren; $i++) {
		my $child_ref = $children_ref->[$i];
		$child_ref->mergeSelfLeaf();
	}
}

# T_bar edge labels for trivial ccs are replaced by the single site of trivial cc
sub simplifyTrivialCCEdgeLabels {
	my $self = shift;
	my $ccs_ref = shift if @_;

	return unless defined ($ccs_ref);

	if ($self->site() > -1) {  # not the root node
		my $site = $self->site();
		if ( (defined $site->[1]) && (scalar(@{$site->[1]}) == 1) ) { # trivial cc label
			croak "Non-trivial CC's cc-restricted sequence is of length 1. Bug!\n" if (scalar(@{$ccs_ref->[$site->[0]]}) > 1);
			$self->site($ccs_ref->[$site->[0]]->[0]);
		}

		my $allsites_ref = $self->allsites();
		for (my $i = 0; $i < scalar(@{$allsites_ref}); $i++) {
			my $item = $allsites_ref->[$i];
			if ( (defined $item->[1]) && (scalar(@{$item->[1]}) == 1) ) { # trivial cc label
				croak "Non-trivial CC's cc-restricted sequence is of length 1. Bug!\n" if (scalar(@{$ccs_ref->[$item->[0]]}) > 1);
				$allsites_ref->[$i] = $ccs_ref->[$item->[0]]->[0];
			}

		}

	 }


	my $size = @{$self->children()};
	my @nodes = @{$self->children()};
	for (my $i = 0; $i < $size; $i++) {
		my $node = $nodes[$i];
		$node->simplifyTrivialCCEdgeLabels($ccs_ref);
	}
}


###########################################################################
package multiTree;

use Exporter;
@ISA = (Exporter);

BEGIN { import multiTree::Node };

use Carp;

sub new {
	my $type = shift;
	my $class = ref($type) || $type;

        my $root = multiTree::Node->new(-1); #root node does not represent any site mutation
	bless {Root => $root}, $class;
}

sub root {
  my ($self, $newroot) = @_;
  $self->{Root} = $newroot if defined $newroot;
  $self->{Root};
}

sub build {
	my ($self, $curNode, $site, $array_ref) = @_;

	#print "curr node = ", $curNode->site(), "\n";
	#print "array = ", @array, "\n";
	#print "site = ", $site, "\n";

	my $size = scalar(@{$array_ref});
	for (my $i = 0; $i < $size; $i++) {
		if ($array_ref->[$i] == $site) {
			my $child = multiTree::Node->new($i); # zero-based site index
			$curNode->addChild($child);
			$self->build($child, $i+1, $array_ref);
		}
	}
}

