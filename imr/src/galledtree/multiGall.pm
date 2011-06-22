package multiGall;
use Carp;
use multiTree;

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

	$self->{rpoint} = undef;
	$self->{rnode} = undef; # stores a 'Tree::Node' object
	$self->{prefix} = []; # stores 'Tree::Node' objects
	$self->{suffix} = []; # stores 'Tree::Node' objects
	$self->{crossOverPts} = []; # stores the sites at which crossover takes place between prefix and suffix

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

sub rpoint {
	my $self = shift;
	$self->{rpoint} = shift if @_;
	return $self->{rpoint};
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

sub crossOverPts {
	my $self = shift;
	$self->{crossOverPts} = shift if @_;
	return $self->{crossOverPts};
}

# computes the on sites based on prefix, suffix and the cross-over points
sub getOnSites
{
	my $self = shift;
	return unless (@_ == 3);

 	my @lastPrefixOns = @{$_[0]};
	my @lastSuffixOns = @{$_[1]};
	my @crossOverPts = @{$_[2]};

	@lastPrefixOns = sort {$a <=> $b} @lastPrefixOns;
	@lastSuffixOns = sort {$a <=> $b} @lastSuffixOns;
	@{$crossOverPts_ref} = sort {$a <=> $b} @{$crossOverPts_ref};

	my @onsites = ();
	my $left; my $right = 0; my $bPrefix = 0;
	my $prefixIndex = 0; my $suffixIndex = 0;

	foreach my $crossOverPt ( @crossOverPts ) {
		$left = $right;
		$right = $crossOverPt;
		$bPrefix = 1 - $bPrefix;

		if ($bPrefix == 1) { # collect prefix ons in range [left, right)
			while ( ($prefixIndex < scalar(@lastPrefixOns)) && ($lastPrefixOns[$prefixIndex] < $left) ) {
				$prefixIndex++;
			}

			while ( ($prefixIndex < scalar(@lastPrefixOns)) && ($lastPrefixOns[$prefixIndex] >= $left)
									&& ($lastPrefixOns[$prefixIndex] < $right) ) {
				push @onsites, $lastPrefixOns[$prefixIndex];
				$prefixIndex++;
			}
		}
		else {  # collect suffix ons in range [left, right)
			while ( ($suffixIndex < scalar(@lastSuffixOns)) && ($lastSuffixOns[$suffixIndex] < $left) ) {
				$suffixIndex++;
			}
			while ( ($suffixIndex < scalar(@lastSuffixOns)) && ($lastSuffixOns[$suffixIndex] >= $left)
									&& ($lastSuffixOns[$suffixIndex] < $right) ) {
				push @onsites, $lastSuffixOns[$suffixIndex];
				$suffixIndex++;
			}
		}
	}

	$left = $right;
	$bPrefix = 1 - $bPrefix;

	if ($bPrefix == 1) { # collect the remaining prefix ons
		while ( $prefixIndex < scalar(@lastPrefixOns) ) {
			push @onsites, $lastPrefixOns[$prefixIndex] if ($lastPrefixOns[$prefixIndex] >= $left);
			$prefixIndex++;
		}
	}
	else {  # collect the remaining suffix ons
		while ( $suffixIndex < scalar(@lastSuffixOns) ) {
			push @onsites, $lastSuffixOns[$suffixIndex] if ($lastSuffixOns[$suffixIndex] >= $left);
			$suffixIndex++;
		}
	}

	return \@onsites;
}

# get the list of 'on' sites of parent. Add my 'on' site to that list.
sub setOnsites {
	my $self = shift;
	return unless @_ == 6;

	my $num_cols; my $parentOns_ref; my $matrix_ref; my $rootSeqID; my $site_counter_ref; my $leaf_counter_ref;
	($num_cols, $parentOns_ref, $matrix_ref, $rootSeqID, $site_counter_ref, $leaf_counter_ref) = @_;

	# walk down prefix and suffix and find 'on' sites for each node
	for (my $side = 0; $side < 2; $side++) {
		my $nodes_ref = ($side == 0) ? $self->{prefix} : $self->{suffix};
		my $numNodes = @{$nodes_ref};

		my $prevOns_ref = $parentOns_ref;

		for (my $i = 0; $i < $numNodes; $i++) {
			my $node = $nodes_ref->[$i];

			if ( $node->findOnSites($num_cols, $prevOns_ref, $matrix_ref, $rootSeqID, $site_counter_ref, $leaf_counter_ref) == -1 ) {
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

	my $crossOverPts_ref = $self->crossOverPts();
	#print "Prefix = @{$lastPrefixOns_ref}\n";
	#print "Suffix = @{$lastSuffixOns_ref}\n";
	#print "Crossovers = @{$crossOverPts_ref}\n";

	my $onsites_ref = $self->getOnSites($lastPrefixOns_ref, $lastSuffixOns_ref, $crossOverPts_ref);
	#print "On sites = @{$onsites_ref}\n";

	# set the final list of 'on' sites to gall
	$self->onsites($onsites_ref);

	#print "On sites of gall are @{$self->onsites()}\n";

	if ( defined $self->sequence() ) {
		# check if the set of 'on' sites at recombination node matches with the corresponding input sequence
		if ( multiMethods->checkSequenceEquality($matrix_ref->[$rootSeqID], $self->onsites(), $matrix_ref->[$self->sequence()]) == 0 ) {
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


	if ( defined $self->sequence() ) { # increment the counter for this leaf sequence
		$leaf_counter_ref->[$self->sequence()]++;
	}

	if ( $self->setOnsites($num_cols, $parentOns_ref, $matrix_ref, $rootSeqID, $site_counter_ref, $leaf_counter_ref) == -1 ) {
		return -1;
	}

	foreach my $child ( @{$self->children()} ) {
		if ( $child->findOnSites($num_cols, \@{$self->onsites()}, $matrix_ref, $rootSeqID, $site_counter_ref, $leaf_counter_ref) == -1 ) {
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
	print $I, "Gall rpoint = ", $self->rpoint() + 1, ", recombination node = ";
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


# returns 1 if the phylogenetic network has a gall with multiple crossovers; 0 otherwise
sub hasMultipleCrossovers {
	my $self = shift;

	my $crossOverPts_ref = $self->crossOverPts();
	return 1 if ( scalar(@{$crossOverPts_ref}) > 1 );

	my $prefixNodes_ref = $self->prefix();
	foreach my $node ( @{$prefixNodes_ref} ) {
		return 1 if ( $node->hasMultipleCrossovers() == 1 );
	}

	my $suffixNodes_ref = $self->suffix();
	foreach my $node ( @{$suffixNodes_ref} ) {
		return 1 if ( $node->hasMultipleCrossovers() == 1 );
	}

	return 0;
}


# every edge label is of the form ( a, b[] ).
sub showSubtreeCompoundEdgeLabels {
	my $self = shift;
	my $indent = shift || 0;
	my $I = ' ' x $indent;
	print $I, "Gall cross-over points = ";

	my $crossOverPts_ref = $self->crossOverPts();
	print $crossOverPts_ref->[0] + 1;
	for (my $i = 1; $i < scalar(@{$crossOverPts_ref}); $i++) {
		print ",", $crossOverPts_ref->[$i] + 1;
	}

	print " recombination node = ";

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
	return $self->rpoint(); # returning some site associated with the gall
				# added to generalize the findGallMissingSitesSequences()
				# method for recombination node and other gall nodes
}


1;
