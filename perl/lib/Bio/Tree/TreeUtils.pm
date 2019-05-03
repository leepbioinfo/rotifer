# $Id: bayesms.pm,v 1.11 2010/07/14 12:46:52 rfsouza Exp $
#
# BioPerl module for Bio::Tree::TreeUtils
#
# Cared for by Robson Francisco de Souza <rfsouza-at-gmail-dot-com>
#
# Copyright Robson Francisco de Souza
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::TreeUtils - A set of subroutines to manipulate and perform
                       calculations on Bio::Tree objects

=head1 SYNOPSIS

  use Bio::TreeIO;
  use Bio::Tree::TreeUtils;
  my $in = Bio::TreeIO->new(-file => 't/data/cat_tre.tre', 
                            -format=>'bayesms');
  my $util = Bio::Tree::TreeUtils->new();
  while( my $tree = $in->next_tree ) {
    $util->remove_nodes($t, 'inparalogs');
  }

=head1 DESCRIPTION

This module implements some subroutines to perform oerations and
calculate statistics from Bio::Tree objects.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Robson Francisco de Souza

Email rfsouza-at-gmail-dot-com

Describe contact details here

=head1 CONTRIBUTORS

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tree::TreeUtils;

use Bio::Root::Root;
use Exporter;
use strict;
use warnings;
use base qw(Exporter Bio::Root::Root);

# Exported methods/tags
our @EXPORT = ();
our @EXPORT_OK = qw();
our %EXPORT_TAGS = ();

=head2 new

 Title   : new
 Usage   : $tu = Bio::Tree::TreeUtils->new()
 Function: create a Bio::Tree::TreeUtils object
 Returns : Bio::Tree::TreeUtils
 Args    : 

=cut

=head2 NODE PROPERTIES

=head2 clean_leaf_names

 Title   : clean_leaf_names
 Usage   : clean_leaf_names($tree, $regexp)
 Function: remove all strings matching $regexp from leaf names
 Returns : Bio::Tree::TreeI
 Args    : Bio::Tree::TreeI and an array of Perl regular expressions

=cut

sub clean_leaf_names {
    my ($self, $tree, @patt) = @_;
    # Changing leaf names
    foreach my $node ($tree->get_nodes) {
	if ($node->is_Leaf) {
	    my $id = $node->id;
	    foreach my $patt (@patt) {
		$id =~ s/$patt//;
		$node->id($id);
	    }
	}
    }
    return $tree;
}

=head2 

 Title   : rename_nodes
 Usage   : rename_nodes($tree, $file)
 Function: Rename nodes in a tree
 Returns : 
 Args    : Bio::Tree::TreeI and
           text table with old names in the 1st column
           and new names on the 2nd column

=cut

sub rename_nodes {
    my $self = shift;
    my $tree = shift;
    my $trans = ref($_[0]) ? $_[0] : _parse_rename_table(@_);

    # Changing leaf names
    foreach my $node ($tree->get_nodes) {
	my $id = $node->id;
	if (defined $id && exists $trans->{$id}) {
	    $node->id($trans->{$id});
	}
    }

    return $tree;
}

sub _parse_rename_table {
    my $file = shift;

    my $hash = {};
    open(TABLE,"<$file") || die "Could not open table $file";
    while (<TABLE>) {
	chomp;
	my @line = split(/\t/);
	$hash->{$line[0]} = $line[1];
    }
    close(TABLE);

    return $hash;
}

=head2 

 Title   : map_node_labels
 Usage   : map_node_labels($tree, $reference, $format, @list)
 Function: Assign nodes in a tree to equivalente nodes in another tree
           The nodes in the target tree are annotated with a 'taxon'
           tag that contains the name of corresponding node in the 
           reference tree. A tag called "leaf_taxa" is also added
           and contains a concatenated list of the taxa for the
           leaves inferred to be the descendents of each node.
 Returns : 
 Args    : Bio::Tree::TreeI object
           reference tree (either a tree file or a Bio::Tree::TreeI)
           reference tree format (undef or Bio::TreeIO parser name)
           (optional) array of best node names from the reference tree

=cut

sub map_node_labels {
    my ($self, $tree, $tax, $format, @targets) = @_;

    # Parse list
    my %targets = ();
    for (my $i=0; $i<=$#targets; $i++) {
	if ( -r $targets[$i] ) {
	    open(TARGET,"<$targets[$i]") || die "Could not open file $targets[$i]";
	    map { chomp; s/^\s*//; s/\s*$//; $targets{$_} = 1 } <TARGET>;
	    close(TARGET);
	} else {
	    $targets{$targets[$i]} = 1;
	}
    }

    # Loading taxonomy
    if ( -r $tax ) {
	$format = 'newick' unless (defined $format);
	my $io = Bio::TreeIO->new(-file=>"<$tax", -format=>$format);
	$tax = $io->next_tree;
	$io->close;
	die "Empty reference tree found while trying to map labels to input tree" unless ($tax->number_nodes > 1);
	$self->postorder_traversal($tax, \&_add_missing_taxon_ids, undef);
    }
    die "map_node_labels needs a tree file or a Bio::Tree::TreeI object as third (function) or second (object) argument."
	unless (ref $tax && $tax->isa("Bio::Tree::TreeI"));

    # Initialize leaf taxonomic annotation
    my @order = $self->postorder_traversal($tree, \&_label2taxonomy, undef, $tax, %targets);

    # Cache old branch properties to avoid losses
    my %support = (); my %length  = ();
    foreach my $node (@order) {
	$support{$node->internal_id} = $node->bootstrap     if (defined $node->bootstrap);
	$length{$node->internal_id}  = $node->branch_length if (defined $node->branch_length);
	$node->branch_length(1); # Use topology
    }

    # Search for better taxon assignments
    my $root = $tree->get_root_node;
    foreach my $node (@order) {
	next if ($node->is_Leaf);
	my $taxon = _get_taxon($tax, $node); # Taxon assigment using original root
	my %taxa  = map { (_get_taxon($tax, $_)->id,1) } grep { $_->is_Leaf } $node->get_all_Descendents;
	next if (!grep { $taxon->id ne $_ } keys %taxa); # No need to revise this node if all descendents are from the same taxon

	# Change root to find the smallest set of leves descending from this node
	foreach my $desc ($node->each_Descendent) {
	    $tree->reroot($desc);
	    my @descs = grep { $_->is_Leaf } $node->get_all_Descendents;
	    my @taxa  = map { _get_taxon($tax, $_) } @descs;
	    if (scalar @taxa < 2) {
		warn "WARNING: apparent loss of data while mapping taxonomy onto tree (".join(" ",scalar(@taxa),map {$_->id} @descs).")\n";
		next if (scalar @taxa < 1);
	    }
	    my ($newTax) = scalar(@taxa) > 1 ? $tax->get_lca(@taxa) : $taxa[0]->ancestor;
	    if ($newTax->height < $taxon->height) { # Lower node in the taxonomy is better!
		%taxa = map { ($_->id,1) } @taxa;
		if (scalar(keys %targets)) {
		    my @selected = $tax->get_lineage_nodes($newTax);
		    my ($selected) = grep { exists $targets{$_->id} } (@selected, $newTax);
		    $newTax = $selected if (defined $selected);
		}
		$taxon = $newTax;
		#print join("\t","Better",$node->id || $node->internal_id, $taxon->id || $taxon->internal_id),"\n";
	    }
	}

	my @taxa = sort keys %taxa;
	$node->set_tag_value('taxon',$taxon->id);
	$node->set_tag_value("leaf_taxa", join(":",@taxa));
	$tree->reroot($root) if ($tree->get_root_node ne $root); # Reset original root
    } # foreach my $node (@order)

    # Reset support values and length: lost by Bioperl!!!!!!
    foreach my $node (@order) {
	my $id = $node->internal_id;
	$node->bootstrap($support{$id})    if (exists $support{$id});
	$node->branch_length($length{$id}) if (exists  $length{$id});
    }
    $root->branch_length('0.0');

    return @order;
}

# Since all leaves must have names, post-order traversal ensures 
# all descendents will have IDs during traversal and this subroutine
# will ensure all nodes in the taxonomy have IDs
sub _add_missing_taxon_ids {
    my ($node, $tree) = @_;

    $node->branch_length(1); # Topology only
    unless (defined $node->id && length $node->id) {
	my @ids = map { $_->id } $node->each_Descendent;
	$node->id( join("|",@ids) );
    }

    return $node;
}

# This subroutine depends on the leaves having names in the format <taxon_id>|<seq_id>
sub _label2taxonomy {
    my ($node, $tree, $taxonomy, %targets) = @_;

    my $tax = undef;
    if ($node->is_Leaf) {
	($tax = $node->id) =~ s/\|.+$//;
	$tax = $taxonomy->find_node(-id => $tax);
	die "Unable to proceed: I cannot find the taxon for ".$node->id." in your taxonomy.\n" unless (defined $tax);
    } else {
	my @taxa = map { $_->get_tag_values('taxon') } $node->each_Descendent;
	@taxa = map { $taxonomy->find_node(-id => $_) } @taxa;
	$tax = scalar @taxa > 1 ? $taxonomy->get_lca(@taxa) : scalar @taxa == 1 ? $taxa[0] : return $node;
    }

    if (scalar(keys %targets)) {
	my @selected = $taxonomy->get_lineage_nodes($tax);
	my ($selected) = grep { exists $targets{$_->id} } (@selected, $tax);
	$tax = $selected if (defined $selected);
    }

    if (defined $tax) {
	$node->set_tag_value('taxon',$tax->id);
	$node->set_tag_value('leaf_taxa',$tax->id);
    }

    return $node;
}

# Retriece Bio::Tree::Node object corresponding to other node's taxon
sub _get_taxon {
    my ($tree, $node) = @_;
    my ($taxID) = $node->get_tag_values('taxon');
    my ($taxon) = $tree->find_node(-id => $taxID);
    ($taxon) = $tree->find_node(-internal_id => $taxID) if (!defined $taxon);
    return $taxon;
}

=head2 TOPOLOGY EXPLORATION

=head2 

 Title   : postorder_traversal
 Usage   : 
 Function: Execute callback subroutine on each node
           during a postorder traversal of the tree

           A postorder traversal is when the deepest 
           subtrees are visited first, i.e. nodes are
           visited in this order:

              +-+ 0
            +-2
            | +-+ 1
           -6
            | +-+ 3
            +-5
              +-+ 4

 Returns : an array of whatever the callback subroutine 
           returns for each node
 Args    : Bio::Tree::TreeI
           (optional) subroutine reference. If undefined,
                      the subroutine will just return the
                      nodes in the order they were visited.
           (optional) first node in subtree (root if undef)
           (optional) other arguments for the call back subroutine

=cut

sub postorder_traversal {
    my ($self, $tree, $subref, $node, @args) = @_;
    $node = $tree->get_root_node if (!defined $node);
    $subref = \&_default_callback unless (defined $subref);

    my @results = ();
    foreach my $desc ($node->each_Descendent) {
	push(@results, $self->postorder_traversal($tree, $subref, $desc, @args));
    }
    push(@results, $subref->($node, $tree, @args));

    return @results;
}

sub _default_callback {
    return defined $_[0] ? $_[0] : $_[1]->get_root_node;
}

=head2 TREE EDITION/SELECTION

=head2 

 Title   : subtree
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub subtree {
    my $self = shift;
    my $tree = shift;
    my @nodeID = @_;

    # Find nodes
    my @nodes = map { $tree->findnode_by_id($_) } @nodeID;
    for (my $i=0; $i<=$#nodeID; $i++) {
	die "Unknown node $nodeID[$i]"
	    unless (ref($nodes[$i]) && $nodes[$i]->isa("Bio::Tree::NodeI"));
    }

    # Reroot and find Last Common Ancestor
    $tree->reroot(pop(@nodes));
    my $lca = $tree->get_lca(-nodes => \@nodes);
    my $anc = $lca->ancestor;
    $tree->reroot($lca);
    $anc->remove_all_Descendents;
   $lca->remove_Descendent($anc);
#    my ($first) = $lca->each_Descendent;
#    $tree->reroot($first);
#    $first->remove_Descendent($lca);

    return $tree;
}

=head2 

 Title   : remove_nodes
 Usage   : remove_nodes($tree, [ @args ])
 Function: 
 Returns : 
 Args    : Bio::Tree::TreeI object
           and reference to an array of 
           node identifiers

           or

           Bio::Tree::TreeI object,
           subroutine name or reference 
           and subroutine arguments

 Supported subroutines are

 inparalogs -> remove all inparalogs (no arguments)

 sametaxon  -> map all nodes to taxons in a reference
               tree and recursively remove leaves so
               as to leave no representatives of the 
               same taxon as direct siblings. For this
               method to work, remove_nodes will be
               called with the syntax

               remove_nodes($tree,'sametaxon',@other)

               where @other are the same arguments as
               those required by L<map_node_labels()>
               AFTER (and excluding) the tree object.

 bootstrap  -> collapse all branches with e-value below
               a threshold, thus generating polytomies
               in the tree

 Notes    : when using 'inparalogs' or 'sametaxon' all
            taxon identifiers should be the first part
            of leaf names, i.e. all leaves should have
            names like

            <taxon>|<sequence identifier>

            Additionally, both 'inparalogs' and 'sametaxon'
            preservwe only the leaves with the smallest
            branch length among siblings of the same taxon.

=cut

sub remove_nodes {
    my ($self, $tree, $array) = @_;

    # Reduce three redundancy based on another tree (taxonomy)
    if ($array->[0] eq 'sametaxon') {
	$self->map_node_labels($tree, @{$array}[1..$#{$array}]);
	return $self->remove_equivalent_branches($tree, \&_sametaxon);
    } 

    # Remove unsupported branches
    elsif ($array->[0] eq 'bootstrap') {
	my @unsupported = ();
	my @keep = ();
	my $root = 0;
	foreach my $node ($tree->get_nodes) {
	    my $b = $node->bootstrap;
	    if (defined $b && length $b && $b < $array->[1]) {
		push(@unsupported,$node->internal_id);
		$root = 1 if ($node eq $tree->get_root_node);
	    } else {
		push(@keep, $node);
	    }
	}
	if ($root) {
	    @keep = sort { $b->is_Leaf <=> $a->is_Leaf } @keep; # Internal nodes have precedence
	    $tree->reroot($keep[0]);
	}
	$tree->splice(-remove_internal_id => \@unsupported, -preserve_lengths=>1) if (scalar @unsupported);
	return @unsupported;
    }

    my @removed = ();
    for (my $i=0; $i<=$#{$array}; $i++) {
	# Compile PERL code once
	if ($array->[$i] =~ /\s*sub\s*\{.+\}\s*/) {
	    $array->[$i] = eval "$array->[$i]";
	    die "Error while compiling subroutine to remove nodes:\n$@" if ($@);
	}

	# Apply filters
	if ($array->[$i] eq 'inparalogs') {
	    push(@removed, $self->remove_equivalent_branches($tree, \&_inparalogs));
	} elsif (ref($array->[$i]) eq "CODE") {
	    push(@removed, $array->[$i]->($tree));
	} else {
	    $tree->splice(-remove_id => [ $array->[$i] ], -preserve_lengths => 1);
	    $self->contract_linear_paths_but_preserve_length($tree);
	    my @internal = grep { !$_->is_Leaf && $_->ancestor } $tree->get_nodes;
	    while (my $in = shift(@internal)) {
		$tree->reroot($in);
		$self->contract_linear_paths_but_preserve_length($tree);
	    }
	}
    }

    return @removed;
}

=head2 

 Title   : remove_equivalent_branches
 Usage   : remove_equivalent_branches($tree, $coderef)
 Function: remove all redundant sister leaves and internal nodes

           Redundancy of sister branches is defined by the subroutine
           $coderef, that will be passed node pairs from the tree and
           should return TRUE whenever the first node is considered
           equivalent to the second.

           The best nodes are the ones whose distance to the last common
           ancestor of a clade composed only by equivalent nodes is the
           shortest (i.e. shorter total branch length to the LCA).

 Returns : number of nodes removed
 Args    : Bio::Tree::TreeI and a reference to a subroutine

=cut

sub remove_equivalent_branches {
    my ($self, $tree, $subref, $preserve) = @_;
    $preserve = 1 unless (defined $preserve);

    my %root = ();
    my %removed = ();
    my $last = -1;
    while ($last != scalar(keys %removed)) {
	$last = scalar(keys %removed);

	# For a given rooted tree, compare all pairs of sibling leafs and maybe remove some (if $subref returns TRUE)
	my @leafs = sort { $a->branch_length <=> $b->branch_length } $tree->get_leaf_nodes;
	foreach my $leaf (@leafs) {
	    next if (exists $removed{$leaf->internal_id});
	    my $ancestor = $leaf->ancestor;
	    my @siblings = sort { $a->branch_length <=> $b->branch_length || $a->id cmp $b->id } grep { $_->is_Leaf } $ancestor->each_Descendent;
	    #print join("\t","LEAF",$leaf->id, $last, scalar(keys %removed), scalar(@leafs), scalar(@siblings), map { $_->id} $ancestor->each_Descendent),"\n";
	    next unless (scalar @siblings); # No leafs that are siblings
 	    for (my $i=0; $i<$#siblings; $i++) {
		next if (exists $removed{$siblings[$i]->internal_id});
		for (my $j=$i+1; $j<=$#siblings; $j++) {
		    next unless ($subref->($siblings[$i], $siblings[$j]));
		    $ancestor->remove_Descendent($siblings[$j]);
		    $removed{$siblings[$j]->internal_id} = 1;
		    #print STDERR join("\t","REMOVED",$siblings[$j]->id,$last,scalar(keys %removed),map { $_->id} $ancestor->each_Descendent),"\n";
		}
	    }
	}

	# Remove single child nodes
	if ($preserve) {
	    $self->contract_linear_paths_but_preserve_length($tree);
	} else {
	    $tree->contract_linear_paths;
	}

	# Nothing was removed: check whether a leaf is a direct descendent of the root and try a new root
	if ($last == scalar(keys %removed)) {
	    #my @internal = grep { !$_->is_Leaf && $_->ancestor && !exists($root{$_->internal_id}) } $tree->get_nodes;
	    #next unless (scalar @internal);
	    #my $newroot = shift(@internal);
	    #$tree->reroot($newroot);
	    #print "NEW ROOT: ",$newroot->internal_id,"\t",scalar(@internal),"\n";
	    #$root{$newroot->internal_id} = 1;
	    #$last = 0;

	    my @internal = $tree->get_root_node->each_Descendent;
	    my @leaves = grep { $_->is_Leaf } @internal;
	    if (scalar @leaves) {
		@internal = grep { !$_->is_Leaf && !exists($root{$_->internal_id}) } @internal;
		next unless (scalar @internal);
		my $newroot = shift(@internal);
		my $ancestor = $newroot->ancestor;
		my $support  = $newroot->bootstrap if (defined $ancestor);
		$tree->reroot($newroot);
		$ancestor->bootstrap($support) if (defined $support && defined $ancestor);
		#print STDERR "NEW ROOT: ",$newroot->id,"\t",scalar(@internal),"\n";
		if ($preserve) {
		    $self->contract_linear_paths_but_preserve_length($tree);
		} else {
		    $tree->contract_linear_paths;
		}
		$root{$newroot->internal_id} = 1;
		$last = 0;
	    }
	} else {
	    %root = ();
	}
    }

    return $last;
}

=head2 

 Title   : contract_linear_paths_but_preserve_length
 Usage   : contract_linear_paths_but_preserve_length($tree, $remove)
 Function: Version of Bio::Tree::TreeFunctionsI::contract_linear_paths
           that preserves branch lengths
 Returns : 
 Args    : Bio::Tree::TreeI and 
           boolean (controls whether the tree should be reroot after
                    contraction)

=cut

sub contract_linear_paths_but_preserve_length {
    my ($self, $tree, $reroot) = @_;
    my @remove = @_;
    foreach my $node ($tree->get_nodes) {
	if ($node->ancestor && $node->each_Descendent == 1) {
	    push(@remove, $node->internal_id);
	}
    }
    $tree->splice(-remove_internal_id => \@remove, -preserve_lengths => 1) if @remove;
    if ($reroot) {
        my $root = $tree->get_root_node;
        my @descs = $root->each_Descendent;
        if (@descs == 1) {
            my $new_root = shift(@descs);
            $tree->set_root_node($new_root);
            $new_root->ancestor(undef);
        }
    }
    return 1;
}

=head2 METRICS

=head2

 Title   : distances
 Usage   : distances($t, $normalize)
 Function: Calculate distances between nodes in a tree
 Returns : two references and one array

           - a reference to a matrix of distances derived by
             summing branch lengths in the path connecting 
             every two nodes

           - a reference to a matrix of distances derived by
             counting the number of internal nodes in the 
             path connecting every two nodes

           - an array of Bio::Tree::Node including all nodes 
             added to the distance matrices (default: 
             all leaves)

 Args    : a Bio::Tree::Node object
           (optional) boolean, whether to normalize distances
           (optional) array of Bio::Tree::Node objects in the
                      tree that you to calculate distances for
                      Default: all leaves (no internal nodes)

=cut

sub distances {
    my ($self, $tree, $normalize, @leaves) = @_;
    @leaves = map { $_->id($_->id || $_->internal_id); $_ } grep { $_->is_Leaf } $tree->get_nodes unless (scalar @leaves);

    # Calculate length
    my $length = $tree->total_branch_length;
    my $nofInternalNodes = $tree->number_nodes - scalar(@leaves);
    #print join(" ",$file,$length),"\n";

    # Calculate distances
    my $dp = [];
    my $nn = [];
    $dp->[0][0] = 0; # Distance to itself is zero
    $nn->[0][0] = 0; # Distance to itself is zero
    for (my $k=1; $k<=$#leaves; $k++) {
	# Loading path from $leaves[$k] to root
	my (@ancsk, %kIDs);
	my $depth = 0;
	my $node  = $leaves[$k];
	while (defined $node) {
	    $kIDs{$node->internal_id} = $depth++;
	    push(@ancsk,$node);
	    $node = $node->ancestor;
	}
	#print join("\t",@kIDs),"\n";

	$dp->[$k][$k] = 0; # Distance to itself is zero
	$nn->[$k][$k] = 0; # Distance to itself is zero
	for (my $l=0; $l<$k; $l++) {
	    my $root  = undef;
	    my @ancs  = @ancsk;
	    my @ancsl = ();
	    my $node  = $leaves[$l];
	    while (defined $node) {
		if (exists $kIDs{ $node->internal_id }) {
		    @ancs = @ancs[0..$kIDs{$node->internal_id}];
		    $root = $node;
		    last;
		}
		unshift(@ancsl, $node);
		$node = $node->ancestor;
	    }
	    @ancs = (@ancs, @ancsl) if (scalar @ancsl);

	    # Add all branch lengths and count nodes in the path connecting k and l
	    my $nofNodes = 0;
	    my $distance = 0;
	    foreach my $node (@ancs) {
		$distance += $node->branch_length unless ($node->internal_id == $root->internal_id);
		$nofNodes++ unless ($node->internal_id == $leaves[$k]->internal_id || $node->internal_id == $leaves[$l]->internal_id);
	    }

#	      print join(" ",
#			 (map { $_->id || $_->internal_id } ($leaves[$k], $leaves[$l], $root)), ":", 
#			 $nofNodes, $distance, $tree->distance(-nodes => [ $leaves[$k], $leaves[$l] ]), ":", 
#			 (map { $_->id || $_->internal_id } @ancs), ":",
#			 (map { $_->id || $_->internal_id } @ancsk), ":",
#			 (map { $_->id || $_->internal_id } @ancsl)
#		  ),"\n";

	    $dp->[$k][$l] = $distance;
	    $dp->[$l][$k] = $distance;
	    $nn->[$k][$l] = $nofNodes;
	    $nn->[$l][$k] = $nofNodes;
	    if ($normalize) {
		$dp->[$k][$l] /= $length;
		$dp->[$l][$k] /= $length;
		$nn->[$k][$l] /= $nofInternalNodes;
		$nn->[$l][$k] /= $nofInternalNodes;
	    }
	} # for (my $l=0; $l<$k; $l++)

    } # for (my $k=1; $k<=$#leaves; $k++)

    return ($dp, $nn, @leaves);
}

=head2 INTERNAL METHODS

=head2 

 Title   : _inparalogs
 Usage   : _inparalogs($node1, $node2)
 Function: Check whether two nodes came from the same organism
           based on the simple naming rule: <taxon>|<leaf name>
 Returns : boolean
 Args    : two Bio::Tree::Node objects or their identifiers

=cut

sub _inparalogs {
    my ($n1,$n2) = @_; 
    $n1 = $n1->id || $n1;
    $n2 = $n2->id || $n2;
    $n1 =~ s/\|\S+//;
    $n2 =~ s/\|\S+//;
    return $n1 eq $n2;
}

=head2 

 Title   : _inparalogs
 Usage   : _inparalogs($node1, $node2)
 Function: Check whether two nodes came from the same organism
           based on the simple naming rule: <taxon>|<leaf name>
 Returns : boolean
 Args    : two Bio::Tree::Node objects or their identifiers

=cut

sub _sametaxon {
    my ($n1,$n2) = @_;
    ($n1) = $n1->get_tag_values('taxon') || '_unknown1';
    ($n2) = $n2->get_tag_values('taxon') || '_unknown2';
    return $n1 eq $n2;
}

1;
