jQuery(document).ready(function(){
  $( ".migration" ).hover(
  function() {
    $( this ).append( $( "<span>What?</span>" ) );
  }, function() {
    $( this ).find( "span:last" ).remove();
  }
);
});
