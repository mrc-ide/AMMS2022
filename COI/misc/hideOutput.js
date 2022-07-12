$(document).ready(function() {
  // Select all <pre> tags that do not have class 'r'
  $output = $("pre").not(".r");
  // Add the show/hide-button to each output chunk
  $output.prepend("<div class=\"showopt\">Show Results</div><br/>");
  // Select the <pre> tags, then choose their <code> child tags and toggle visibility
  $output.children("code").css({display: "none"});

  // When the show/hide-button is clicked, toggle the current state and
  // change the button text
  $(".showopt").click(function() {
    $btn = $(this);
    $chunk = $(this).parent().children("code");
    if($btn.html() === "Show Results") {
      $btn.html("Hide Results");
    } else {
      $btn.html("Show Results");
    }
    $chunk.slideToggle("fast", "swing");
  });
});
