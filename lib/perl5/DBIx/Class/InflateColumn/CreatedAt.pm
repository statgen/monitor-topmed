package DBIx::Class::InflateColumn::CreatedAt;

use base qw(DBIx::Class);

use DateTime;
use DateTime::Format::MySQL;

sub insert {
  my $self = shift;
  my $rec  = $self->next::method(@_);
  my $now  = DateTime->now();

  $rec->update(
    {
      created_at => DateTime::Format::MySQL->format_datetime($now),
    }
  );

  return $rec;
}

1;
