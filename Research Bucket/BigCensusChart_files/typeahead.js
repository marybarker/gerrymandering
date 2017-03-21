(function($) {
$(document).ready(function(){
    $( "#searchbox, #mobile_searchbox" ).catcomplete({
        source: function( request, response ) {
            $.ajax({
                url: "/suggest", // Needs to set up apache forwarding to hit SOA typeahead url or JSP/Java forwarding
                dataType: "xml",
                data: {
                    query: request.term,
                    stateGeo: $("stateGeo").val(),
                    operationName: "httpsearch"
                },
                success: function( data ) {
                    // Merge the two types we want to display and iterate over them
                    response( $.map($.merge($(data).find("Answer"), $(data).find("result")), function(item) {
                        // Only the answer xml has a value component
                        if($(item).find("Value").text() != ""){ // for the answers before results
                            return {
                                ansValue: $(item).find("Value").text(),
                                ansGeography: $(item).find("Geography").text(),                                                                
                                prefix: $(item).find("Prefix").text(),
                                suffix: $(item).find("Suffix").text(),
                                ansLastUpdated: $(item).find("LastUpdated").text(),
                                // Don't grab the system desc
                                ansDescription: $(item).find("Description").first().text(),
                                ansSystemDesc: $(item).find("System").find("Description").text(),
                                ansSystemURL: $(item).find("System").find("URL").text()
                            }
                        }else if($(item).attr("name") != ""){ // for regular typeahead results
                            return {
                                label: $(item).attr("name"),
                                value: $(item).attr("name"),
                                category: $(item).attr("category")
                            }
                        }
                    }));
                    
                }
            });
        },
        minLength: 3,
        select: function( event, ui ) {
            if(ui.item){
                $('#searchbox, #mobile_searchbox').val(ui.item.value);
                $('input[name=cssp]').val('Typeahead');
            }
            
            $('#searchbox, #mobile_searchbox').closest("form").submit();
        },
        messages: {
            noResults: '',
                results: function() {}
        },
        open: function() {
//          $(this).dialog("widget").css("z-index", 10000);
            $( this ).removeClass( "ui-corner-all" ).addClass( "ui-corner-top" );
        },
        close: function() {
            $( this ).removeClass( "ui-corner-top" ).addClass( "ui-corner-all" );
        }
    }); 
});

// Custom Typeahead
$.widget( "custom.catcomplete", $.ui.autocomplete, {
    _renderMenu: function( ul, items ) {
        var self = this, currentCategory = "";
        
        $.each( items, function( index, item ) {
            // double check for blank nodes from the jquery find method
            if(item.ansValue != null || item.label != null){
                // Decorate the system field to include latest estimate date
                var sysName = "";
                
                if(item.ansLastUpdated != ""){
                    var d = new Date(item.ansLastUpdated * 1000); // Given seconds originally
                    sysName = d.getFullYear() + " " + item.ansSystemDesc;
                }else{
                    sysName = item.ansSystemDesc;
                }
                
                if(item.ansValue != null){

                        ul.append('<li class="autocomplete-instant-answer"> ' +
                                '<div class="ia-label"><span class="ia-geo">' + item.ansGeography + '</span>&nbsp|&nbsp<span class="ia-stat">' + item.ansDescription + '</span></div>' +
                                '<div class="ia-stat-value">'+ item.prefix+'' +item.ansValue.replace(/(\d)(?=(\d\d\d)+(?!\d))/g, "$1,") +' ' +item.suffix+'</div>' +
                                
                                '<div class="ia-system"><a href="' + item.ansSystemURL  + '?cssp=Typeahead">Source: '+ sysName + '</a></div>' +
                                '</li>');
                    
                }else if ( item.category != currentCategory ) {
                    ul.append( "<li class='ui-autocomplete-category'>" + item.category + "</li>" );
                    currentCategory = item.category;
                }
                
                self._renderItemData( ul, item );
            }
            
        });
    }
});
})($acn);