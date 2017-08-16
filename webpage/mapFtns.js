function myMap(){
  
  map = new google.maps.Map(document.getElementById("mapCanvas"), {
    zoom:7, 
    center: new google.maps.LatLng(centroid[0],centroid[1]) 
  });
  
  map.data.loadGeoJson(jsonfilename, {}, function(features){yet_to_assign = features.length;});
  
  map.data.forEach(function(feature){
    feature.setProperty('District', 'Not Assigned');
  });
  
  map.data.setStyle(function(feature){
    myColor = allColors[feature.getProperty('District')];
    return ({
      fillColor: myColor, 
      strokeWeight: .1
    });
  });
  
  map.data.addListener('click', function(event){
    // get district this feature was in before click
    var value = event.feature.getProperty("District");

    // set the feature to be in district currentDistrict
    event.feature.setProperty("District", currentDistrict);
    
    // update total unassigned features.
    if (currentDistrict != "Not Assigned"){
      yet_to_assign--;
    }
    if(value >= 0){
      yet_to_assign++;
    }
    document.getElementById("unassigned").innerHTML="Unassigned VTDS: "+yet_to_assign.toString();
    
    // now update stats on current district
    var statsString='';
    for(var name of list_of_functions_to_compute){
      statsString += (calculateAll(currentDistrict, name)) + '<br>';
    }
    document.getElementById("stats").innerHTML = statsString;
  });
  
  map.data.addListener('mouseover', function(event){
    map.data.revertStyle();
    map.data.overrideStyle(event.feature, {strokeWeight:2});
  });
  
  map.data.addListener('mouseout', function(event){
    map.data.revertStyle();
  });
}

