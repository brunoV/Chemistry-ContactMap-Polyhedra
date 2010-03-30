use strict;
use warnings;
use Test::More;
use Chemistry::MacroMol;
use Chemistry::File::PDB;
use Bio::Tools::Run::QCons;

BEGIN { use_ok( 'Chemistry::ContactMap::Polyhedra' ) }

my @mols = map { Chemistry::MacroMol->new } (0..1);

my $cmap = Chemistry::ContactMap::Polyhedra->new( structures => \@mols );

isa_ok $cmap, 'Chemistry::ContactMap::Polyhedra';

my @methods = qw(residue_contacts atom_contacts add_bonds_to_structures);
can_ok( $cmap, @methods);

my $mol1 = Chemistry::MacroMol->read('t/data/A.pdb');
my $mol2 = Chemistry::MacroMol->read('t/data/B.pdb');

$cmap = Chemistry::ContactMap::Polyhedra->new( structures => [ $mol1, $mol2 ] );

is ref $cmap->contacts, 'ARRAY';

isa_ok( $cmap->contacts->[0], 'Chemistry::Bond' );

# compare results against those returned by QCons
my $q = Bio::Tools::Run::QCons->new( file => 't/data/Q.pdb', chains => ['A', 'B'] );

is_deeply( $cmap->residue_contacts, $q->residue_contacts );
# residue contacts don't match because PerlMol reindexes atom numbers

is scalar @{$cmap->contacts}, scalar @{$q->atom_contacts};

done_testing();
