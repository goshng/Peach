package Gall;
use Carp;
use Tree;

sub new {
	my $type = shift;
	my $class = ref($type) || $type;
	my $self = {};

	$self->{parent} = undef;
	$self->{seekParent} = 0; # states whether or not to direct the edge to parent. 0 is undirected, 1 is directed.
	$self->{children} = []; # these lead from the recombination node
				# stores 'Gall' and 'Tree::Node' objects
	$self->{seekChildren} = {};  # states whether or not to direct the edges to each chid. 0 is undirected, 1 is directed.

	$self->{sequence} = undef; # sequence of recombination node

	$self->{onsites} = []; # sites that will be 1 in the recombination node
	$self->{zeroOneState} = []; # stores 011000xxx.. the state of 0s and 1s in node

	$self->{rInterval} = [];
	$self->{rnode} = undef; # stores a 'Tree::Node' object
	$self->{prefix} = []; # stores 'Tree::Node' objects
	$self->{suffix} = []; # stores 'Tree::Node' objects

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

# sequence of recombination node
sub sequence {
	my $self = shift;
	$self->{sequence} = shift if @_;
	return $self->{sequence};
}

sub rInterval {
	my $self = shift;
	$self->{rInterval} = shift if @_;
	return $self->{rInterval};
}

# creates an array of prefix nodes from an array of prefix sites(mere numbers)
sub prefix {
	my $self = shift;
	return $self->{prefix} unless @_;

	my $prefix_ref = $self->{prefix};

	foreach my $site (@_) {
		push @{$prefix_ref}, Tree::Node->new($site);
	}

	return $self->{prefix};
}

# creates an array of suffix nodes from an array of suffix sites(mere numbers)
sub suffix {
	my $self = shift;
	return $self->{suffix} unless @_;

	my $suffix_ref = $self->{suffix};

	foreach my $site (@_) {
		push @{$suffix_ref}, Tree::Node->new($site);
	}

	return $self->{suffix};
}

sub getPrefixSites {
	my $self = shift;
	my @preSites = ();

	my $prefix_ref = $self->{prefix};

	foreach my $node (@{$prefix_ref}) {
		push @preSites, $node->site();
	}

	return @preSites;
}

sub getSuffixSites {
	my $self = shift;
	my @sufSites = ();

	my $suffix_ref = $self->{suffix};

	foreach my $node (@{$suffix_ref}) {
		push @sufSites, $node->site();
	}

	return @sufSites;
}

# returns the first prefix site. returns -1 if no prefix exists
sub getFirstPrefixSite {
	my $self = shift;
	my $prefix_ref = $self->{prefix};

	foreach my $node (@{$prefix_ref}) {
		return $node->site();
	}

	return -1;
}

# returns the first suffix site. returns -1 if no suffix exists
sub getFirstSuffixSite {
	my $self = shift;
	my $suffix_ref = $self->{suffix};

	foreach my $node (@{$suffix_ref}) {
		return $node->site();
	}

	return -1;
}

# adds only one child, even if multiple children are presented
sub addChild {
	my $self = shift;
	if ( @_ > 0 ) {
		my $child = shift;

		#print "Adding...\n";
		#$child->display();
		#print "as child of ...\n";
		#$self->display();

		push @{$self->{children}}, $child;
		$child->parent($self);
	}
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

# get the list of 'on' sites of parent. Add my 'on' site to that list.
sub setOnsites {
	my $self = shift;
	return unless @_ == 8;

	my $num_cols; my $parentOns_ref; my $oldNewSiteMap_ref; my $duplicateColsHash_ref;
	my $matrix_ref; my $rootSeqID; my $site_counter_ref; my $leaf_counter_ref;
	($num_cols, $parentOns_ref, $oldNewSiteMap_ref, $duplicateColsHash_ref, $matrix_ref, $rootSeqID, $site_counter_ref, $leaf_counter_ref) = @_;

	# walk down prefix and suffix and find 'on' sites for each node
	for (my $side = 0; $side < 2; $side++) {
		my $nodes_ref = ($side == 0) ? $self->{prefix} : $self->{suffix};
		my $numNodes = @{$nodes_ref};

		my $prevOns_ref = $parentOns_ref;

		for (my $i = 0; $i < $numNodes; $i++) {
			my $node = $nodes_ref->[$i];

			if ( $node->findOnSites($num_cols, $prevOns_ref, $oldNewSiteMap_ref, $duplicateColsHash_ref,
							$matrix_ref, $rootSeqID, $site_counter_ref, $leaf_counter_ref) == -1 ) {
				return -1;
			}
			$prevOns_ref = \@{$node->onsites()};
		}
	}

	# take the 'on' sites of last prefix and last suffix if they exist
	my $prefixNodes_ref = $self->{prefix};
	my $numPrefixNodes = @{$prefixNodes_ref};
	my $lastPrefixOns_ref = [];

	if ( $numPrefixNodes > 0 ) {
		my $lastPrefixNode = $prefixNodes_ref->[$numPrefixNodes - 1];
		$lastPrefixOns_ref = \@{$lastPrefixNode->onsites()};
	}
	else { # no prefix nodes, take 'on sites' from the coalescent node
		my $parent = $self->parent();
		croak "Gall must have a coalescent node. Bug!\n" if (not defined $parent);
		$lastPrefixOns_ref = \@{$parent->onsites()};
	}

	my $suffixNodes_ref = $self->{suffix};
	my $numSuffixNodes = @{$suffixNodes_ref};
	my $lastSuffixOns_ref = [];

	if ( $numSuffixNodes > 0 ) {
		my $lastSuffixNode = $suffixNodes_ref->[$numSuffixNodes - 1];
		$lastSuffixOns_ref = \@{$lastSuffixNode->onsites()};
	}
	else { # no suffix nodes, take 'on sites' from the coalescent node
		my $parent = $self->parent();
		croak "Gall must have a coalescent node. Bug!\n" if (not defined $parent);
		$lastSuffixOns_ref = \@{$parent->onsites()};
	}

	# remove P_R and S_L
	$self->rInterval()->[1] = $oldNewSiteMap_ref->[$self->rInterval()->[1]];  # set to corresponding old recombination interval
	my $rPoint = $self->rInterval()->[1];

	my @lastPrefixOns = @{$lastPrefixOns_ref};
	my @lastSuffixOns = @{$lastSuffixOns_ref};

	for (my $i = 0; $i < scalar(@lastPrefixOns); $i++) {
		if ( $lastPrefixOns[$i] >= $rPoint ) {
			splice @lastPrefixOns, $i, 1;
			$i--;
		}
	}

	for (my $i = 0; $i < scalar(@lastSuffixOns); $i++) {
		if ( $lastSuffixOns[$i] < $rPoint ) {
			splice @lastSuffixOns, $i, 1;
			$i--;
		}
	}

	# set the final list of 'on' sites to gall
	my @selfOns = (@lastPrefixOns, @lastSuffixOns);

	$self->onsites(\@selfOns);

	#print "On sites of gall interval = [", $self->rInterval()->[0], ",", $self->rInterval()->[1], "] are @{$self->onsites()}\n";

	if ( defined $self->sequence() ) {
		# check if the set of 'on' sites at recombination node matches with the corresponding input sequence
		if ( Methods->checkSequenceEquality($matrix_ref->[$rootSeqID], $self->onsites(), $matrix_ref->[$self->sequence()]) == 0 ) {
			return -1; # if not equal
		}

		$self->zeroOneState($matrix_ref->[$self->sequence()]);
	}
	else {
		my $zeroOneState_ref = Methods->getZeroOneState($matrix_ref->[$rootSeqID], $num_cols, $self->onsites());
		$self->zeroOneState($zeroOneState_ref);
	}
}


sub findOnSites {
	my $self = shift;
	return unless @_ == 8;

	my $num_cols; my $parentOns_ref; my $oldNewSiteMap_ref; my $duplicateColsHash_ref;
	my $matrix_ref; my $rootSeqID; my $site_counter_ref; my $leaf_counter_ref;
	($num_cols, $parentOns_ref, $oldNewSiteMap_ref, $duplicateColsHash_ref, $matrix_ref, $rootSeqID, $site_counter_ref, $leaf_counter_ref) = @_;

	if ( defined $self->sequence() ) { # increment the counter for this leaf sequence
		$leaf_counter_ref->[$self->sequence()]++;
	}

	if ( $self->setOnsites($num_cols, $parentOns_ref, $oldNewSiteMap_ref, $duplicateColsHash_ref,
				$matrix_ref, $rootSeqID, $site_counter_ref, $leaf_counter_ref) == -1 ) {
		return -1;
	}

	foreach my $child ( @{$self->children()} ) {
		if ( $child->findOnSites($num_cols, \@{$self->onsites()}, $oldNewSiteMap_ref, $duplicateColsHash_ref,
						$matrix_ref, $rootSeqID, $site_counter_ref, $leaf_counter_ref) == -1 ) {
			return -1;
		}
	}
}


# leaves correspond to the matrix relabelled with a original leaf as root. So, relabelled again to correspond to the original leaf
sub reLabelLeaves {
	my $self = shift;
	return unless @_ == 1;

	my $leafReorder_ref = $_[0];

	my $zeroOneState_ref = $self->zeroOneState();
	my $numSites = @{$zeroOneState_ref};

	for (my $i = 0; $i < $numSites; $i++) {
		if ($leafReorder_ref->[$i] == 1) {
			$zeroOneState_ref->[$i] = 1 - $zeroOneState_ref->[$i];
		}
	}

	# walk down prefix and suffix and relabel the leaves there
	for (my $side = 0; $side < 2; $side++) {
		my $nodes_ref = ($side == 0) ? $self->{prefix} : $self->{suffix};
		my $numNodes = @{$nodes_ref};

		for (my $i = 0; $i < $numNodes; $i++) {
			my $node = $nodes_ref->[$i];
			$node->reLabelLeaves($leafReorder_ref);
		}
	}

	foreach my $child ( @{$self->children()} ) {
		$child->reLabelLeaves($leafReorder_ref);
	}
}



sub parent {
	my $self = shift;
	$self->{parent} = shift if @_;
	return $self->{parent};
}


# increasing sites and sequence numbers by 1 to be 1-based instead of being traditionally 0-based
sub showSubtree {
	my $self = shift;
	my $indent = shift || 0;
	my $I = ' ' x $indent;
	print $I, "Gall interval = [", $self->rInterval()->[0] + 1, ",", $self->rInterval()->[1] + 1, "] recombination node = ";
	print @{$self->zeroOneState()}, " (row ";
	(defined $self->sequence()) ? print $self->sequence() + 1 : print "?"; # check if sequence exists
	print ")\n";
	print $I, "  Prefix side\n";
	my $prefixNodes_ref = $self->prefix();
	foreach my $node ( @{$prefixNodes_ref} ) {
		$node->showSubtree($indent+2);
	}

	print $I, "  Suffix side\n";
	my $suffixNodes_ref = $self->suffix();
	foreach my $node ( @{$suffixNodes_ref} ) {
		$node->showSubtree($indent+2);
	}

	print $I, "  Below recombination node\n";
	my $childNodes_ref = $self->children();
	foreach my $node ( @{$childNodes_ref} ) {
		$node->showSubtree($indent+2);
	}
}

# every edge label is of the form ( a, b[] ).
sub showSubtreeCompoundEdgeLabels {
	my $self = shift;
	my $indent = shift || 0;
	my $I = ' ' x $indent;
	print $I, "Gall interval = [", $self->rInterval()->[0] + 1, ",", $self->rInterval()->[1] + 1, "] recombination node = ";
	print @{$self->zeroOneState()}, " (row ";
	(defined $self->sequence()) ? print $self->sequence() + 1 : print "?"; # check if sequence exists
	print ")\n";
	print $I, "  Prefix side\n";
	my $prefixNodes_ref = $self->prefix();
	foreach my $node ( @{$prefixNodes_ref} ) {
		$node->showSubtreeCompoundEdgeLabels($indent+2);
	}

	print $I, "  Suffix side\n";
	my $suffixNodes_ref = $self->suffix();
	foreach my $node ( @{$suffixNodes_ref} ) {
		$node->showSubtreeCompoundEdgeLabels($indent+2);
	}

	print $I, "  Below recombination node\n";
	my $childNodes_ref = $self->children();
	foreach my $node ( @{$childNodes_ref} ) {
		$node->showSubtreeCompoundEdgeLabels($indent+2);
	}
}

# returns recombination point (which is also a gall site)
sub site {
	my $self = shift;
	return $self->rInterval()->[1]; # returning some site associated with the gall
				# added to generalize the findGallMissingSitesSequences()
				# method for recombination node and other gall nodes
}


1;
