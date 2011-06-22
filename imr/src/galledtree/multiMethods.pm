package multiMethods;
use Carp;
use multiTree;
use multiGall;


# if possible, constructs a perfect phylogeny and returns pointer to the root of the phylogeny
# else returns -1
sub buildPerfPhylogeny
{
	my $self = shift;
	return unless @_ == 1;

	my $matrix_ref = $_[0];

	my $num_rows = @{$matrix_ref};
	my $num_cols = @{$matrix_ref->[0]};

	#print "Matrix is \n";
	#for (my $row = 0; $row < $num_rows; $row++) {
	#	print "@{$matrix_ref->[$row]}\n";
	#}

	# remove duplicate columns. Add later to get edges with multiple labels
	my $uniqueCol_ref, $duplicateColsHash_ref;
	($uniqueCol_ref, $duplicateColsHash_ref) = $self->getUniqueColsDuplicatesList($matrix_ref);

	my $oldNewSiteMap_ref = $self->getOldNewSiteMap($duplicateColsHash_ref, $num_cols);

	$matrix_ref = $self->getUniqueColumnsMatrix($matrix_ref, $uniqueCol_ref);
	$num_cols = @{$matrix_ref->[0]};

	# sort (descending order) the columns treating a column as a binary number
	# exponentiation operator is **

	my @decimals;
	for (my $col = 0; $col < $num_cols; $col++) {
		my $value = 0;
		for (my $row = 0; $row < $num_rows; $row++) {
			$value += $matrix_ref->[$row][$col] * (2 ** ($num_rows - $row - 1) );
		}

		push @decimals, $value;
	}

	#print "decimals = @decimals\n";

	my @siteOrder = (0 .. $num_cols-1);
	@siteOrder = sort { $decimals[$b] <=> $decimals[$a] } @siteOrder;
	#print "Site order = @siteOrder\n";

	# construct matrix M' using input matrix M by reordering the columns as in @order
	my @ordMatrix;
	for (my $row = 0; $row < $num_rows; $row++) {
		for (my $col = 0; $col < $num_cols; $col++) {
   			$ordMatrix[$row][$col] = $matrix_ref->[$row][$siteOrder[$col]];
		}
	}

	#print "Ordered matrix is\n";
	#for (my $row = 0; $row < $num_rows; $row++) {
	#	print "@{$ordMatrix[$row]}\n";
	#}

	# construct matrix L(i,j) using the matrix M'
	my @modMatrix;
	for (my $row = 0; $row < $num_rows; $row++) {
		for (my $col = $num_cols - 1; $col >= 0; $col--) {

			if ($ordMatrix[$row][$col] == 0) {
				 $modMatrix[$row][$col] = 0;
			}
			else
			{
				$modMatrix[$row][$col] = 0;
				for (my $i = $col - 1; $i >= 0 ; $i--) {
					if ($ordMatrix[$row][$i] == 1) {
						#print "Row = $row, Col = $col, I = $i\n";
						$modMatrix[$row][$col] = $i + 1;
						last;
					}
				}
			}

		}
	}

	# construct hash C_i using the matrix M'. Used to assign leaf sequences
	# map a site to a sequence index for fast acccess
	my %seqOrder;
	for (my $row = 0; $row < $num_rows; $row++) {
		my $col;
		for ($col = $num_cols - 1; $col >= 0; $col--) {
			if ($ordMatrix[$row][$col] == 1) {
				$seqOrder{$col} = $row;
				last;
			}
		}

		if ($col == -1) { # could not find a '1' in the 'root' row
			$seqOrder{-1} = $row; # hard-coding -1 since the site of 'root' node is -1
		}
	}

	#print "Sequence order ...\n";
	#my @seqKeys = keys (%seqOrder);
	#foreach my $key (@seqKeys) {
	#	print "Site $key => sequence $seqOrder{$key}\n";
	#}


	#print "Modified matrix is\n";
	#for (my $row = 0; $row < $num_rows; $row++) {
	#	print "@{$modMatrix[$row]}\n";
	#}

	# For every column j, set L(j) equal to the largest modMatrix(i,j) such that ordMatrix(i,j) = 1
	my @maxModColumn;
	for (my $col = 0; $col < $num_cols; $col++) {
		my $max = 0;

		for (my $row = 0; $row < $num_rows; $row++) {
			if ($ordMatrix[$row][$col] == 1) {
				$max = $modMatrix[$row][$col] if ($modMatrix[$row][$col] > $max);
			}
		}

		push @maxModColumn, $max;
	}

	#print "Maximum column values L(j) = @maxModColumn\n";

	my $perfect = 1; # flag to determine if perfect phylogeny exists. Assumed to exist initially.

	# verify if a perfect phylogeny is possible, that is if modMatrix(i,j) = maxModColumn(j) for every ordMatrix(i,j) = 1
	for (my $col = 0; ($col < $num_cols) && ($perfect == 1); $col++) {
		for (my $row = 0; ($row < $num_rows) && ($perfect == 1); $row++) {
			if ($ordMatrix[$row][$col] == 1) {
				if ($modMatrix[$row][$col] != $maxModColumn[$col]) {
					#print "Perfect phylogeny is not possible($row, $col)\n";
					$perfect = 0;
				}
			}
		}
	}

	if ($perfect == 1) {
		my $tree = multiTree->new();
		my $root = $tree->root();

		$tree->build($root, 0, \@maxModColumn);
		$root->addSubtreeLeaves(\%seqOrder);
		$root->reLabelSites(\@siteOrder); #maps sites to unique columns matrix

		#$root->showSubtree();
		#print "After deleting sequence less subtrees...\n";

		$root->removeSequencelessSubTree(-1);
		#$root->showSubtree();

		$root->reLabelSites($oldNewSiteMap_ref); # remaps sites to duplicate columns matrix
		$root->addDuplicateEdgeLabels($duplicateColsHash_ref);

		return $root;
	}
	else { # it is possible that sequences below nodes collected by bottom-up traversal do not form a perfect phylogeny.
	       # since sequences rightfully below nodes in other galls were collected here.
		return -1;
	}
}


# accepts an array reference and returns the reference with duplicates removed and sorted in non-descending order
sub removeDuplicates {
	my $self = shift;
	return unless @_ == 1;

	my $list_ref = shift;
	@list = sort {$a <=> $b} @{$list_ref};
	$list_ref = \@list;

	for (my $i = 0; $i < scalar(@{$list_ref}) - 1; $i++) {
		if ( $list_ref->[$i] == $list_ref->[$i+1] ) { # adjacent identical entries => duplicates
			splice @{$list_ref}, $i, 1;
			$i--;
		}
	}

	return $list_ref;
}


sub getZeroOneState {
	my $self = shift;
	return unless @_ == 3;

	my $zeroOneState_ref = [];

	my $rootSeq_ref; my $num_cols; my $onSites_ref;
	($rootSeq_ref, $num_cols, $onSites_ref) = @_;

	my @onSites = sort {$a <=> $b} @{$onSites_ref};

	for (my $i = 0; $i < $num_cols; $i++) {
		my $state = $rootSeq_ref->[$i];

		if ( scalar(@onSites) > 0  ) {
			if ( $onSites[0] == $i) {
				splice @onSites, 0, 1;
				$state = 1 - $state;
			}
		}

		push @{$zeroOneState_ref}, $state;
	}

	return $zeroOneState_ref;
}


# compares the list of 'on' sites with the alleged corresponding input sequence for equality
# Note: We are not checking if every site mutates only once in the network. However, if a site
# mutates twice and the corresponding nodes have a parent-child relationship, then the duplicates
# also get propagated, resulting in a fail at this method.
sub checkSequenceEquality {
	my $self = shift;
	return unless @_ == 3;

	my $rootSeq_ref; my $onSites_ref; my $sequence_ref;
	($rootSeq_ref, $onSites_ref, $sequence_ref) = @_;

	#print "Checking sequence equality ...\n";
	#print "root sequence = @{$rootSeq_ref}\n";
	#print "onsites = @{$onSites_ref}\n";
	#print "sequence = @{$sequence_ref}\n";

	my @onSites = sort {$a <=> $b} @{$onSites_ref};
	#print "onsites = @onSites\n";

	my $onSitesIndex = 0;

	for (my $sequenceIndex = 0; $sequenceIndex < scalar(@{$sequence_ref}); $sequenceIndex++) {
		if ( $sequence_ref->[$sequenceIndex] != $rootSeq_ref->[$sequenceIndex] ) {
			if ( $onSitesIndex == scalar(@onSites) ) {
				#print "Not equal 3\n";
				return 0;
			}


			if ( $onSites[$onSitesIndex] != $sequenceIndex ) {
				#print "Not equal 1\n";
				return 0;
			}

			$onSitesIndex++;
		}
	}

	if ( $onSitesIndex != scalar(@onSites) ) {
		#print "Not equal 2\n";
		return 0;
	}

	#print "Equal\n";
	return 1;
}


# accepts a reference to an array of counters to leaf occurences in galled tree and checks if
# all leaves occur exactly once; If so, returns 1 (success) else returns 0 (failure)
sub leafExactlyOnce
{
	my $self = shift;
	return unless @_ == 1;

	my $vector_ref = shift;

	for (my $i = 0; $i < scalar(@{$vector_ref}); $i++) {
		if ( $vector_ref->[$i] > 1 ) {
			print "Row ", $i+1, " of unique row matrix occuring more than once.\n";
			return 0;
		}

		if ( $vector_ref->[$i] == 0 ) {
			print "Row ", $i+1, " of unique row matrix missing.\n";
			return 0;
		}
	}

	return 1;
}

# accepts a reference to an array of counters to site occurences in galled tree and checks if
# all sites mutate no more than once; If so, returns 1 (success) else returns 0 (failure)
sub siteAtmostOnce
{
	my $self = shift;
	return unless @_ == 1;

	my $vector_ref = shift;

	for (my $i = 0; $i < scalar(@{$vector_ref}); $i++) {
		if ( $vector_ref->[$i] > 1 ) {
			print "Site ", $i+1, " of unique column matrix mutating more than once.\n";
			return 0;
		}
	}

	return 1;
}


# accepts a reference to a matrix, a reference to an array of sites in some
# CC and returns a reference to a sub-matrix with only those sites/columns
sub getCCMatrix
{
	my $self = shift;
	return unless @_ == 2;

	my $matrix_ref; my $cc_ref;;
	($matrix_ref, $cc_ref) = @_;

	my $num_rows = @{$matrix_ref};
	my $num_cols =  @{$cc_ref};

	my $subMatrix_ref;

	foreach my $j ( @{$cc_ref} ) {
		for (my $i = 0; $i < $num_rows; $i++) {
			push @{$subMatrix_ref->[$i]}, $matrix_ref->[$i][$j];
		}
	}

	return $subMatrix_ref;
}


# returns a hash of the unique row values(as string). Each hash value is the reference to an array of row numbers having that unique row value.
sub groupIdenticalRows
{
	my $self = shift;
	return unless @_ == 1;

	my $input_matrix_ref = $_[0];

	my %uniqueRows;
	my $num_rows = @{$input_matrix_ref};
	my $num_cols = @{$input_matrix_ref->[0]};
	my %hashRows;

	for (my $row = 0; $row < $num_rows; $row++) {
		my $newtext = "";
		for (my $col = 0; $col < $num_cols; $col++) {
			$newtext .= $input_matrix_ref->[$row][$col];
		}

		push @{$uniqueRows{$newtext}}, $row;
	}

	#my $num_unique_rows = scalar( keys(%uniqueRows) );
	#print "The unique $num_unique_rows rows and their corresponding row numbers are\n";
	#foreach my $item (sort { $a <=> $b } (keys % uniqueRows) ) {
	#	print "$item,--> @{$uniqueRows{$item}} \n";
	#}
	#print "\n";

	return \%uniqueRows;
}


# adds a new column to given matrix. The column has a 1 in the rows specified and a 0 in the rest.
sub addColumn
{
	my $self = shift;
	return unless @_ == 3;

	my $matrix_ref; my $rows_ref; my $num_rows;
	($matrix_ref, $rows_ref, $num_rows) = @_;

	# assuming that the specified row numbers are sorted
	my $rowID = 0; my $numRowIDs = @{$rows_ref};
	for (my $i = 0; $i < $num_rows; $i++) {
		if ( ($rowID < $numRowIDs)  && ($rows_ref->[$rowID] == $i) ) {
			push @{$matrix_ref->[$i]}, 1;
			$rowID++;
		}
		else {
			push @{$matrix_ref->[$i]}, 0;
		}
	}
}


# returns two hash references. The first refers to a hash with unique columns as keys
# The second refers to a hash with keys as only those columns with duplicates. Maps key to an array of duplicate columns to its right.
sub getUniqueColsDuplicatesList
{
	my $self = shift;
	return unless @_ == 1;

	my $input_matrix_ref = $_[0];

	my %hashCols;

	my %uniqueCols; # keys are unique columns (leftmost chosen when duplicates exist)
	my %duplicateCols; # keys are only those unique columns for which duplicates exist. Maps a key to an array of its duplicate columns
	my $num_rows = @{$input_matrix_ref};
	my $num_cols = @{$input_matrix_ref->[0]};

	for (my $col = 0; $col < $num_cols; $col++) {
		my $newtext = "";
		for (my $row = 0; $row < $num_rows; $row++) {
			$newtext .= $input_matrix_ref->[$row][$col];
		}

		if ( defined $hashCols{$newtext} ) {
			push @{$duplicateCols{$hashCols{$newtext}}}, $col;
		}
		else {
			$hashCols{$newtext} = $col; # leftmost column hashes to a column value
			$uniqueCols{$col} = 'c'; # store unique columns as keys in a hash
		}
	}

	return (\%uniqueCols, \%duplicateCols);
}


sub getUniqueColumnsMatrix
{
	my $self = shift;
	return unless @_ == 2;

	my $original_matrix_ref = $_[0];
	my %uniqueCols = %{$_[1]};
	my $num_rows = @{$original_matrix_ref};

	my @uniqueMatrix;

	for (my $rowItem = 0; $rowItem < $num_rows; $rowItem++ ) {
		my $newtext = "";

		foreach my $colItem (sort { $a <=> $b } (keys % uniqueCols) ) {
			$newtext .= $original_matrix_ref->[$rowItem][$colItem];
			$newtext .= " ";
		}

		my (@row) = split (/\s+/, $newtext);
		push @uniqueMatrix, \@row;
	}

	return \@uniqueMatrix;
}


# returns an array of the indices of unique rows in matrix
sub getUniqueRows
{
	my $self = shift;
	return unless @_ == 1;

	my $input_matrix_ref = $_[0];

	my %uniqueRows;
	my $num_rows = @{$input_matrix_ref};
	my $num_cols = @{$input_matrix_ref->[0]};
	my %hashRows;

	for (my $row = 0; $row < $num_rows; $row++) {
		my $newtext = "";
		for (my $col = 0; $col < $num_cols; $col++) {
			$newtext .= $input_matrix_ref->[$row][$col];
		}

		if ( ! $hashRows{$newtext} ) {
			#push @uniqueRows, $row;
			$uniqueRows{$row} = 'c'; # store unique rows as keys in a hash instead of array (for faster look ups)
		}

		$hashRows{$newtext} = 'c'; # some key-value pair
	}

	#my $num_unique_rows = scalar( keys(%uniqueRows) );
	#print "The unique $num_unique_rows rows are ";
	#foreach my $item (sort { $a <=> $b } (keys % uniqueRows) ) {
	#	print "$item ";
	#}
	#print "\n";

	return \%uniqueRows;
}


sub getUniqueRowsMatrix
{
	my $self = shift;
	return unless @_ == 2;

	my $original_matrix_ref = $_[0];
	my %uniqueRows = %{$_[1]};

	my @uniqueMatrix;

	foreach my $rowItem (sort { $a <=> $b } (keys % uniqueRows) ) {
		push @uniqueMatrix, $original_matrix_ref->[$rowItem];
	}

	#my $num_rows = @uniqueMatrix;
	# printing each row of the matrix
        #for (my $i = 0; $i < $num_rows; $i++) {
	#	print "@{$uniqueMatrix[$i]}\n";
	#}
	#print "\n";

	return \@uniqueMatrix;
}


# returns a map (array) of sites in unique columns matrix to those in unique rows matrix
sub getOldNewSiteMap
{
	my $self = shift;
	return unless @_ == 2;

	my $duplicateColsHash_ref; my $num_cols;
	($duplicateColsHash_ref, $num_cols) = @_;

	my @oldNewSiteMap = (0..$num_cols-1);

	for (my $i = 0; $i < scalar(@oldNewSiteMap); $i++) {
		if ( defined $duplicateColsHash_ref->{$oldNewSiteMap[$i]} ) {
			foreach my $j ( @{$duplicateColsHash_ref->{$oldNewSiteMap[$i]}}  ) {
				for (my $k = $i + 1; $k < scalar(@oldNewSiteMap); $k++) {
					if ( $oldNewSiteMap[$k] == $j ) {
						splice @oldNewSiteMap, $k, 1;
						last;
					}
				}
			}
		}
	}

	return \@oldNewSiteMap;
}


sub recomputeCCconflicts
{
	my $self = shift;
	return unless @_ == 4;

	my $confColHash_ref; my $ccs_ref; my $recombPoints_ref; my $oldNewSiteMap_ref;
	($confColHash_ref, $ccs_ref, $recombPoints_ref, $oldNewSiteMap_ref) = @_;

	my %hash;

	my $numDistinctSites = @{$oldNewSiteMap_ref};

	for (my $i = 0; $i < $numDistinctSites; $i++) {
		$hash{$oldNewSiteMap_ref->[$i]} = $i;
	}

	my @ccs; my @recombPoints;

	my $numOldCCs = @{$ccs_ref};
	for (my $i = 0; $i < $numOldCCs; $i++) {
		my @row;

		foreach my $j (@{$ccs_ref->[$i]}) {
			push @row, $hash{$j} if defined $hash{$j};
		}

		if ( scalar(@row) > 0 ) {
			push @ccs, \@row;
			if ( scalar(@row) == 1 ) {
				push @recombPoints, -1;
			}
			else {
				if ( not defined $hash{$recombPoints_ref->[$i]} ) {
					croak "Leftmost suffix site not being retained in unique columns matrix. Wrong!!!\n";
				}
				push @recombPoints, $hash{$recombPoints_ref->[$i]};
			}
		}
	}

	#print "New ccs are...\n";
	#foreach my $aref (@ccs) {
	#	print "@{$aref}\n";
	#}

	my %confColHash;
	my $num_cols = scalar(keys %{$confColHash_ref});
	for (my $i = 0; $i < $num_cols; $i++) {
		if ( defined $hash{$i} ) {
			my $ihash = $hash{$i};
			$confColHash{$ihash} = [];

			foreach my $j ( @{$confColHash_ref->{$i}} ) {
				push @{$confColHash{$ihash}}, $hash{$j} if defined $hash{$j};
			}
		}
	}

	#print "New conflicts are...\n";
	#$num_cols = scalar(keys % confColHash);
	#for (my $i = 0; $i < $num_cols; $i++) {
	#	print "Column $i conflicts with the columns @{$confColHash{$i}}\n";
	#}
	#print "\n";

	return (\%confColHash, \@ccs, \@recombPoints);
}

# returns a shallow copy of matrix
sub getShallowMatrixCopy {
	my $self = shift;
	return unless (@_ == 1);

	my $given_ref = $_[0];
	my $rowCount = @{$given_ref};

	my @copyMatrix;
	for (my $i = 0; $i < $rowCount; $i++) {
		push @copyMatrix, $given_ref->[$i];
	}

	return \@copyMatrix;
}

# returns a deep copy of matrix
sub getDeepMatrixCopy {
	my $self = shift;
	return unless (@_ == 1);

	my $given_ref = $_[0];
	my $return_ref;

	my $rowCount = @{$given_ref};
	my $colCount = @{$given_ref->[0]};

	for (my $i = 0; $i < $rowCount; $i++) {
		my @dupRow;

		for (my $j = 0; $j < $colCount; $j++) {
			push @dupRow, ${$given_ref->[$i]}[$j];
		}


		push @{$return_ref}, \@dupRow;
	}

	return $return_ref;
}

1;
