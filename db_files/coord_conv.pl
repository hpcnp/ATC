#!/usr/bin/perl
##
#
use Astro::Coords;
use Astro::Telescope;
use Astro::Catalog;
use Time::Piece;

my $time = new Time::Piece;


# adjust for time zone??
# $time += 6*3600;

$c = new Astro::Coords( planet => 'jupiter' );

# Location
# 9.956318, -83.995503, 1480m
$c->telescope( new Astro::Telescope( Name => 'Home', Long => -117.72303, Lat => 33.6974747, Alt => 22 ));

$c->datetime( $time );

$az = $c->az( format => 'd' );
$el = $c->el( format => 'd' );
$ra = $c->ra( format => 'd' );
$dec = $c->dec( format => 'd' );
print "az: $az\nel: $el\n\nra: $ra\ndec: $dec\n";
