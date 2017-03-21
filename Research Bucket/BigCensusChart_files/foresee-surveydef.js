(function(){

    var FSR;

    // Do we support AMD?
    var supports_amd =
        typeof(window.define) === 'function' && window.define.amd &&
            (!window.FSR || window.FSR.supportsAMD);

    if(!supports_amd)
        FSR = window.FSR;
    else
        FSR = {};

    FSR.surveydefs = [{
    name: 'browse',
    section: 'low-volume',
    invite: {
        when: 'onentry'
    },
    pop: {
        when: 'later'
    },
    criteria: {
        sp: 100,
        lf: 2
    },
    include: {
        urls: ['econ(?!(((.)*www(.)*))|((.)*census02(.)))', 'census.gov/govs/', 'census.gov/manufacturing/', 'census.gov/construction/', 'census.gov/const/', 'census.gov/retail/', 'census.gov/wholesale/', 'census.gov/services/', '2010.census.gov/partners/', 'census.gov/schools/2010_census/', 'census.gov/mtis/', 'census.gov/cgi-bin/briefroom/BriefRm', 'census.gov/ces/']
    }
}, {
    name: 'browse',
    section: 'high-volume',
    invite: {
        when: 'onentry'
    },
    pop: {
        when: 'later'
    },
    criteria: {
        sp: 75,
        lf: 2
    },
    include: {
        urls: ['.']
    }
}]
FSR.properties = {
    repeatdays: 0,

    repeatoverride : false,

    altcookie : {
    },

    language : {
        locale : 'en'
    },

    exclude : {
    
        variables: [{
            name: 'fsr$ip',
            value: [/148\.129\.*\.*/, /172\.16\.0\.*/, /192\.168\.0\.*/]
        }]
    },
	/* Invite branding sample property
    brands : [{"c":"Foresee","p":33}, {"c":"Answers", "p":33}, {"c":"ForeseeByAnswers", "p":33}],
	*/
    zIndexPopup : 10000,

    ignoreWindowTopCheck : false,

    ipexclude : 'fsr$ip',

    mobileHeartbeat : {
        delay : 60, /*mobile on exit heartbeat delay seconds*/
        max : 3600  /*mobile on exit heartbeat max run time seconds*/
    },

    invite : {

        // For no site logo, comment this line:
        siteLogo : "sitelogo.gif",

        //alt text fore site logo img
		siteLogoAlt : "",

        /* Desktop */
        dialogs : [[{
            reverseButtons: false,
            headline: "We'd welcome your feedback!",
            blurb: "Thank you for visiting the US Census Bureau. You have been selected to participate in a brief customer satisfaction survey to let us know how we can improve your experience.",
            noticeAboutSurvey: "The survey is designed to measure your entire experience, please look for it at the <u>conclusion</u> of your visit.",
            attribution: "This survey is conducted by an independent company ForeSee, on behalf of the site you are visiting.",
            closeInviteButtonText: "Click to close.",
            declineButton: "No, thanks",
            acceptButton: "Yes, I'll give feedback",
            error: "Error",
            warnLaunch: "this will launch a new window"

        }]],

        exclude : {
            urls:[],
            referrers:[],
            userAgents:[],
            browsers:[],
            cookies:[],
            variables:[]
			// [name (content), http-equiv (content), itemprop (content),  charset] possible attributes for meta tag element http://devdocs.io/html/meta
            // metas:[{"name":{"key":"value", "content":"value"}}, {"http-equiv":{"key":"value", "content":"value"}}, {"itemprop":{"key":"value", "content":"value"}}, {"charset":{"key":"value"}}]
        
        },
        include : {
            local : [ '.' ]
        },

        delay : 0,
        timeout : 0,

        hideOnClick : false,

        hideCloseButton : false,

        css : 'foresee-dhtml.css',

        hide : [],

        hideFlash: false,

        type : 'dhtml',
        /* desktop */
        // url: 'invite.html'
        /* mobile */
        url : 'invite-mobile.html',
        back: 'url'

        //SurveyMutex: 'SurveyMutex'
    },

    tracker : {
        width : '690',
        height : '415',
        timeout : 3,
		//pu: false,
        adjust : true,
        alert : {
            enabled : true,
            message : 'The survey is now available.'
        },
        url : 'tracker.html'
    },

    survey : {
        width : 690,
        height : 600
    },

    qualifier : {
        footer : '<div id=\"fsrcontainer\"><div style=\"float:left;width:80%;font-size:8pt;text-align:left;line-height:12px;\">This survey is conducted by an independent company ForeSee,<br>on behalf of the site you are visiting.</div><div style=\"float:right;font-size:8pt;\"><a target="_blank" title="Validate TRUSTe privacy certification" href="//privacy-policy.truste.com/click-with-confidence/ctv/en/www.foreseeresults.com/seal_m"><img border=\"0\" src=\"{%baseHref%}truste.png\" alt=\"Validate TRUSTe Privacy Certification\"></a></div></div>',
        width : '690',
        height : '500',
        bgcolor : '#333',
        opacity : 0.7,
        x : 'center',
        y : 'center',
        delay : 0,
        buttons : {
            accept : 'Continue'
        },
        hideOnClick : false,
        css : 'foresee-dhtml.css',
        url : 'qualifying.html'
    },

    cancel : {
        url : 'cancel.html',
        width : '690',
        height : '400'
    },

    pop : {
        what : 'survey',
        after : 'leaving-site',
        pu : false,
        tracker : true
    },

    meta : {
        referrer : true,
        terms : true,
        ref_url : true,
        url : true,
        url_params : false,
        user_agent : false,
        entry : false,
        entry_params : false,
		viewport_size: false,
        document_size: false,
        scroll_from_top: false,
		invite_URL: false
    },

    events : {
        enabled : true,
        id : true,
        codes : {
            purchase : 800,
            items : 801,
            dollars : 802,
            followup : 803,
            information : 804,
            content : 805
        },
        pd : 7,
        custom : {}
    },

    previous : false,

	analytics : {
		google_local : false,
		google_remote : false
	},

    cpps: {
        Econ: {
            source: 'url',
            patterns: [{
                regex: '/econ/index/|/econ/geo/|/econ/survey/|/econ/accomodation/|econ/concentration|econ/economywide/|/econ/other/|/econ/progoverview/|/econ/construction/|/econ/manufacturing/|/econ/retail/|econ/services/|econ/wholesale|/econ/www/|/econ/industry|/econ/isp',
                value: '1'
            }, {
                regex: '/econ/census/|/econ/census02/|/econ/census07/',
                value: '2'
            }, {
                regex: '/econ/cbp/|/econ/nonemployer/',
                value: '4'
            }, {
                regex: '/econ/sbo/|/econ/qfr|econ/susb/|/econ/aces|econ/ict/|/csd/',
                value: '5'
            }, {
                regex: '/epcd/',
                value: '6'
            }, {
                regex: '/foreign-trade/',
                value: '7'
            }, {
                regex: '/govs/',
                value: '8'
            }, {
                regex: '/manufacturing/|/mcd/',
                value: '9'
            }, {
                regex: '/retail/|/wholesale/|/services/|/mtis/|/svsd/',
                value: '10'
            }, {
                regex: '/naics/naicsrch|/eos/',
                value: '11'
            }, {
                regex: '/economic-indicators/',
                value: '12'
            }, {
                regex: '/ces/',
                value: '13'
            }, {
                regex: '/construction/|/const/',
                value: '14'
            }, {
                regex: '.',
                value: '99'
            }]
        },
        delivery_src: {
            source: 'url',
            init: 'none',
            patterns: [{
                regex: 'source=govdelivery',
                value: 'govdelivery'
            }]
        }
    },

    mode : 'first-party'
};

    if(supports_amd)
        define(function(){ return FSR; })
})();
