/* Unlink.js - Remove links to the current page.
 * Mon Mar 26 13:17:02 EDT 2007
 * Lou Scoras <louis DOT j DOT scoras AT census DOT gov>
 * Requires 'JQuery'
 *
 * This script will unlink anchor elements that refer to the current page in
 * the browser.  The canonical use for this pertains to navigation menus.
 * Normally when the user clicks on a menu element, the selected page should
 * remain in the element, however it should not be selectable.
 *
 * This script will remove any such links within the scope of a given ancestor
 * element.  Both the links and the container elements are specified via css
 * selectors.  Furthermore, an optional callback function can be provided
 * which will reverence the container element.  This is useful for doing
 * further styling operations.
 *
 * For instance, say you have a menu which contains links within an unordered
 * list.
 *
 * + ul { class => menu } 
 *   |
 *   +-- li
 *   |   |
 *   |   a { href => foo }
 *   |
 *   +-- li
 *   |   |
 *   |   a { href => current_page }
 *   |
 *   +-- li
 *       |
 *       a { href => bar }
 *   
 * If we want to remove the self references from the menu, a new unlinker
 * object should be created like so:
 *
 *   menuUnlinker = new Unlinker('ul.menu', 'a')
 *
 * The first parameter tells the unlinker how to find the container element
 * for the links.   These are the things we want to style.  Specifying 'a' for
 * the second parameter says to search for all links below it to be candidates
 * for removal.
 *
 * One of the useful things about the unlinker object is that it allows you to
 * specify a callback for styling purposes.  Say whenever links are removed in
 * the above structure, the li element containing the link should be styled
 * differently.  We could specify it like so,
 *
 *   menuUnlinker = new Unlinker('ul.menu > li', 'a', function(li) {
 *     li.addClass('selected');
 *   })g
 *
 * where addClass is a jquery function which addes the 'selected' CSS class to
 * the li tag.
 *
 */


/* The Unlinker object takes two css selectors: the first is for the
 * containing element; the second tells the code how to find the link elements
 * to remove.
 *
 * The callback parameter is optional.  It will be passed the container
 * element whenever a link is removed.
 */
function Unlinker(pq, lq, callback) {
  this.parentQuery = pq;
  this.linkQuery   = lq;

  if (callback) {
    this.callback = callback;
  }
  else {
    this.callback = function(p) {return void(0);};
  }
}

Unlinker.prototype = {
 
  /* Removes a single link from the page.  The link text will be reinserted
   * into the page as a span tag
   */
  removeLink : function(link, parentElement) {
    var matcher = new URLMatcher(link.href);
    if (matcher.isThisPage()) {
      var replacement = document.createElement("span");
      replacement.appendChild(document.createTextNode(link.childNodes[0].nodeValue));
      link.parentNode.replaceChild(replacement, link);
      this.callback(parentElement);
    }
  },

  /* Call this function to kick off the the link removal process
   */
  removeAllLinks : function() {
    var groups = $(this.parentQuery);
    for (var i = 0; i < groups.length; i++) {
      this.removeLinkGroup(groups[i]);
    }
  },

  /* Preforms the processing for each link container element as specified in
   * the constructor
   */
  removeLinkGroup : function(parentElement) {
    var links = $(this.linkQuery, parentElement);
    for (var i = 0; i < links.length; i++) {
      this.removeLink(links[i], parentElement);
    }
  }
};

/* The URLMatcher is used by the unlink code to determine wheter two links are
 * actually equal for this particular bit of functionality.  Relative links
 * are converted into full urls, and any anchor information is ignored.
 */
function URLMatcher(url) {
  this.url = url;
}

URLMatcher.prototype = {

  /* Returns whether or not the matcher considers this link to be the same as
   * the current page
   */
  isThisPage : function() {
    return (document.location == this.canonicalURL());
  },

  /* Trys to normalize the url
   */
  canonicalURL : function() {
    var url = this.url.replace(/#.*$/, '');
    if (this.isRelative()) {
      var pathname = document.location.pathname;
      var filePart = pathname.replace(new RegExp("/[^/]*$"),('/' + url));
  
      return ("http://" + document.location.host + filePart);
    }
    else return url;
  },

  /* Decides whether the url is absolute or relative
   */
  isRelative : function() {
    if (this.url.indexOf(':') > 0) {
      return false;
    }
    else if (this.url.charAt(0) == '/') {
      return false
    }
    return true;
  }
};
