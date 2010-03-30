package Chemistry::ContactMap::Polyhedra;

# ABSTRACT: Calculate contacts with the Voroni Polyhedra method

use Mouse;
use MouseX::Types::Mouse qw(Bool);
use Bio::Tools::Run::QCons;
use Chemistry::MacroMol;
use Chemistry::Bond;
use File::Temp;

with 'Chemistry::ContactMap';

has add_bonds_to_structures => ( is => 'ro', isa => Bool, default => 0 );

has _qcons => (
    is         => 'ro',
    lazy_build => 1,
    handles => [qw( residue_contacts atom_contacts )],
);

sub _build__qcons {
    my $self = shift;

    # Logic:
    # get input structures and write a temporary file with both.
    # Since Qcons needs the molecules to have distinct chain names, we
    # clone the structure, assign arbitrary chain names, and use those
    # to run Qcons with the two modified structures written into a
    # temporary file.

    my ($mol1, $mol2) = map { $_->clone } @{$self->structures};

    foreach my $residue ($mol1->domains) {
        $residue->attr('pdb/chain_id', 'A');
    }

    foreach my $residue ($mol2->domains) {
        $residue->attr('pdb/chain_id', 'B');
    }

    my $temp_structure = Chemistry::MacroMol->new();
    $temp_structure->add_domain($_) for ($mol1->domains, $mol2->domains);

    my $tmp_file     = File::Temp->new;
    my $tmp_filename = $tmp_file->filename . '.pdb';

    $temp_structure->write($tmp_filename);

    my $qcons = Bio::Tools::Run::QCons->new(
        file   => $tmp_filename,
        chains => ['A', 'B'],
    );

    # trigger calculation here: if $tmp_file runs out of scope, it's
    # going to be deleted.
    $qcons->residue_contacts;

    return $qcons;
}

sub _build_contacts {

    # Create Chemistry::Bonds between the atoms that are
    # interacting non covalently. Extra properties are stored in the
    # generic 'attr' attribute. Types are: V, H, I and S for
    # van der Waals, Hydrogen Bond, Ionic pair and Salt bridge,
    # respectively. Properties are dependent upon type, but all of them
    # have an 'area' attribute.

    my $self = shift;

    my ( $mol1, $mol2 ) = @{ $self->structures };

    # Associating atoms from the QCons output and the working
    # MacroMol structure turned out to be quite trickier than I first
    # thought. The only way to make an unambiguous matching is through
    # their residue's 'pdb/sequence_number' attribute. So I first build
    # an array indexed by this attribute. Then, when I want to get the
    # atoms that are participating in the contact, I get the parent
    # residue and fetch the atom by its name (atom names are unique within
    # each residue, I checked it).

    my ( @residues1, @residues2 );

    $residues1[ $_->attr('pdb/sequence_number') ] = $_ for $mol1->domains;
    $residues2[ $_->attr('pdb/sequence_number') ] = $_ for $mol2->domains;

    my @contacts;

    foreach my $contact ( @{ $self->atom_contacts } ) {

        my ($atom1) =
          grep { $_->name eq $contact->{atom1}->{name} }
          $residues1[ $contact->{atom1}->{res_number} ]->atoms;

        my ($atom2) =
          grep { $_->name eq $contact->{atom2}->{name} }
          $residues2[ $contact->{atom2}->{res_number} ]->atoms;

        my $bond = Chemistry::Bond->new(
            type  => $$contact{type},
            order => 0,
            atoms => [ $atom1, $atom2 ],
            attr  => { area => $$contact{area}, },
        );

        foreach my $attr ( keys %$contact ) {
            next if $attr ~~ /atom|type|area/;    # Already added.
            $bond->attr( $attr, $contact->{$attr} );
        }

        push @contacts, $bond;

        if ( $self->add_bonds_to_structures ) {
            $_->add_bond($bond) for ( $mol1, $mol2 );
        }
    }

    return \@contacts;
}

1;

__END__

=head1 SYNOPSIS

    my $cmap = Chemistry::ContactMap::Polyhedra->new( structures => [ $mol1, $mol2 ] );

    foreach my $contact ( @{ $cmap->contacts } ) {

        # get the intervening atoms
        my ($atom1, $atom2) = $contact->atoms;

        say $_->name, $_->number for ($atom1, $atom2);

        say $contact->length(), $contact->attr('area'), $contact->type();
    }

=head1 DESCRIPTION

This class implements the calculation of inter-atom and inter-residue
non-covalent contacts between macromolecules using the Voronoi Polyhedra
Method.

It consumes the L<Chemistry::ContactMap> interface, interacts seamlessly
with L<PerlMol>, and uses L<Bio::Tools::Run::QCons>, a Perl wrapper over
Tsai's QCons program.

=attr structures

Two L<Chemistry::MacroMol> objects whose contact information you want to
calculate. Required.

=attr add_bonds_to_structures

A boolean value that sets whether to add the atom-atom contacts found
between the macromolecules for each atom pair, as a Chemistry::Bond
object of order one. It is set as 0 by default.

=attr contacts

An array reference of L<Chemistry::Bond> objects that represent the
non-covalent (order 0) bonds between atom-atom pairs that are found to
be interacting.

The contact type (Van der Waals, Ionic, Hydrogen bond or Hydrophobic),
the total contact area and the bond length are saved in the object.

See L<Chemistry::Bond> for details on how to use the bond class.

=attr residue_contacts

An array reference with residue-residue contact information. This data
structure is identical to that returned by L<Bio::Tools::Run::QCons>,
see the docs there to find out how to extract the relevant information
from it.

=attr atom_contacts

An array reference with atom-atom contact information. Identical to that
returned by L<Bio::Tools::Run::QCons>. See its doc for more information.
