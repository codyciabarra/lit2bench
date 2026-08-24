/* site.js -- interaction for the Lit2Bench site.
 *
 * No framework and no bundler: this is served straight off GitHub Pages next to
 * index.html. Everything here is progressive -- the page reads and works with
 * JavaScript switched off, and each block below only enhances what's already
 * rendered in the HTML.
 */
(function () {
  'use strict';

  /* ------------------------------------------------------------- release
   * The *fallback* download details -- what renders immediately, and what the
   * page keeps if the API call further down never lands. The live answer comes
   * from GitHub's Releases API at page load, so this constant drifting behind a
   * tag is no longer fatal.
   *
   * It still has to name an asset that exists. `latest/download/<asset>` is a
   * GitHub redirect to the newest release's copy of that filename, so it follows
   * new releases on its own -- but the filename carries the version, and a name
   * no release ever published 404s. That is exactly what shipped at v0.1.0: the
   * site advertised Lit2Bench-0.1.0.dmg before any release was tagged, so the
   * button was dead. Bump `version` and `asset` together when you tag.
   */
  var RELEASE = {
    version: '0.2.0',
    size: '~1.6 MB',
    repo: 'https://github.com/codyciabarra/lit2bench',
    asset: 'Lit2Bench-0.2.0.dmg'
  };
  RELEASE.url = RELEASE.repo + '/releases/latest/download/' + RELEASE.asset;

  /* ----------------------------------------------------------- analytics
   * Cloudflare Web Analytics: visitor counts with no cookies, no localStorage
   * and no cross-site identifiers, which is why it needs no consent banner.
   *
   * Paste the beacon token from the Cloudflare dashboard below to switch it on.
   * While it's empty nothing is injected at all -- no script tag, no request, no
   * console error -- so an unconfigured site is simply a site without analytics
   * rather than a site with a broken one.
   *
   * It has to be the manual beacon rather than Cloudflare's automatic setup:
   * automatic injection only works for proxied (orange-cloud) hostnames, and
   * these records are deliberately DNS-only so GitHub can hold the TLS cert.
   */
  var ANALYTICS_TOKEN = '';
  if (ANALYTICS_TOKEN) {
    var beacon = document.createElement('script');
    beacon.defer = true;
    beacon.src = 'https://static.cloudflareinsights.com/beacon.min.js';
    beacon.setAttribute('data-cf-beacon', JSON.stringify({ token: ANALYTICS_TOKEN }));
    document.head.appendChild(beacon);
  }

  var reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
  var $  = function (s, r) { return (r || document).querySelector(s); };
  var $$ = function (s, r) { return Array.prototype.slice.call((r || document).querySelectorAll(s)); };

  /* --------------------------------------------------------------- download */
  var dlBtn = $('#dl-btn');
  if (dlBtn) {
    dlBtn.href = RELEASE.url;
    $('#dl-version').textContent = RELEASE.version;
    $('#dl-size').textContent = RELEASE.size;

    // Visitors on Windows or Linux shouldn't be handed a .dmg without warning:
    // the app runs fine from source on their machine, only the installer is
    // macOS-only so far.
    var ua = navigator.userAgent;
    if (!/Mac/i.test(ua) || /iPhone|iPad/i.test(ua)) {
      $('#dl-label').textContent = 'Download the macOS .dmg';
      $('#dl-sub').textContent =
        'macOS only for now — on Windows or Linux, run from source (see right).';
    }

    /* Ask GitHub what the newest release actually is, and correct the button.
     *
     * The hardcoded RELEASE above is the fallback, not the source of truth. It
     * has to be *a* real asset name because that's what renders before this
     * request lands (and all that renders if it never does) -- but when the
     * constant drifts behind a tag, as it did at v0.1.0 when the site advertised
     * an asset no release had ever published, this repairs the button instead of
     * serving a 404 to every visitor.
     *
     * Deliberately not awaited and never fatal: the download link is already
     * usable when this fires, so a rate-limited or offline API leaves the page
     * exactly as it was. GitHub allows 60 unauthenticated calls an hour per IP,
     * which is per-visitor here, so one call per page load is well inside it.
     */
    fetch('https://api.github.com/repos/codyciabarra/lit2bench/releases/latest',
          { headers: { Accept: 'application/vnd.github+json' } })
      .then(function (r) { return r.ok ? r.json() : null; })
      .then(function (rel) {
        if (!rel || !rel.tag_name) return;
        var dmg = (rel.assets || []).filter(function (a) {
          return /\.dmg$/i.test(a.name);
        })[0];
        if (!dmg) return;   // a release with no .dmg (source-only tag): leave the fallback

        dlBtn.href = dmg.browser_download_url;
        $('#dl-version').textContent = String(rel.tag_name).replace(/^v/, '');
        // Real bytes beat a hand-maintained "~1.5 MB" that nobody remembers to bump.
        $('#dl-size').textContent = '~' + (dmg.size / 1048576).toFixed(1) + ' MB';
      })
      .catch(function () { /* offline or rate-limited -- fallback already rendered */ });
  }
  var year = $('#year');
  if (year) year.textContent = new Date().getFullYear();

  /* ------------------------------------------------------------------- nav */
  var nav = $('#nav');
  var onScroll = function () { nav.classList.toggle('scrolled', window.scrollY > 8); };
  onScroll();
  addEventListener('scroll', onScroll, { passive: true });

  /* ----------------------------------------------------------------- theme */
  var themeBtn = $('#theme');
  if (themeBtn) {
    themeBtn.addEventListener('click', function () {
      var next = document.documentElement.dataset.theme === 'light' ? 'dark' : 'light';
      document.documentElement.dataset.theme = next;
      try { localStorage.setItem('l2b-theme', next); } catch (e) {}
    });
  }

  /* ---------------------------------------------------------------- reveal */
  var revealables = $$('.reveal');
  if (reduced || !('IntersectionObserver' in window)) {
    revealables.forEach(function (el) { el.classList.add('in'); });
  } else {
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (!e.isIntersecting) return;
        // Stagger siblings so a row of cards arrives as a wave, not a slab.
        var peers = $$('.reveal', e.target.parentNode);
        e.target.style.transitionDelay = (Math.min(peers.indexOf(e.target), 6) * 60) + 'ms';
        e.target.classList.add('in');
        io.unobserve(e.target);
      });
    }, { rootMargin: '0px 0px -12% 0px', threshold: 0.08 });
    revealables.forEach(function (el) { io.observe(el); });

    // Belt and braces: some renderers (headless Chrome without a compositor,
    // and a few in-app webviews) never deliver intersection callbacks. Reveal
    // everything unconditionally after a moment so a missed callback can't
    // leave a section permanently invisible.
    setTimeout(function () {
      revealables.forEach(function (el) { el.classList.add('in'); });
    }, 2500);
  }

  /* -------------------------------------------------------------- counters */
  // Numbers count up once, the first time they scroll into view. The element's
  // text is already correct in the HTML, so this is decoration only.
  $$('[data-count]').forEach(function (el) {
    var target = parseInt(el.dataset.count, 10);
    if (reduced || isNaN(target) || target === 0 || !('IntersectionObserver' in window)) return;
    var obs = new IntersectionObserver(function (entries) {
      if (!entries[0].isIntersecting) return;
      obs.disconnect();
      var t0 = performance.now(), dur = 900;
      (function tick(now) {
        var p = Math.min((now - t0) / dur, 1);
        el.textContent = Math.round(target * (1 - Math.pow(1 - p, 3)));
        if (p < 1) requestAnimationFrame(tick);
      })(t0);
    }, { threshold: 0.6 });
    obs.observe(el);
  });

  /* ---------------------------------------------------------- tool filters */
  var grid = $('#tool-grid');
  if (grid) {
    var tools = $$('.tool', grid);
    var chips = $$('.chip');
    var search = $('#tool-search');
    var empty = $('#no-results');
    var group = 'all';

    var apply = function () {
      var q = (search.value || '').trim().toLowerCase();
      var shown = 0;
      tools.forEach(function (t) {
        var okGroup = group === 'all' || t.dataset.group === group;
        var okText = !q || t.textContent.toLowerCase().indexOf(q) !== -1;
        var show = okGroup && okText;
        t.classList.toggle('hidden', !show);
        if (show) shown++;
      });
      empty.hidden = shown > 0;
    };

    chips.forEach(function (chip) {
      chip.addEventListener('click', function () {
        group = chip.dataset.filter;
        chips.forEach(function (c) { c.setAttribute('aria-pressed', String(c === chip)); });
        apply();
      });
    });
    search.addEventListener('input', apply);

    /* Pointer spotlight: feed the cursor position to the card's ::before
       gradient. Coordinates are per-card and only updated while hovering, so
       this costs nothing when the pointer is elsewhere. */
    if (!reduced && matchMedia('(hover: hover)').matches) {
      tools.forEach(function (t) {
        t.addEventListener('pointermove', function (e) {
          var r = t.getBoundingClientRect();
          t.style.setProperty('--mx', (e.clientX - r.left) + 'px');
          t.style.setProperty('--my', (e.clientY - r.top) + 'px');
        });
      });
    }
  }

  /* ------------------------------------------------------------------ copy */
  $$('.copy').forEach(function (btn) {
    btn.addEventListener('click', function () {
      var code = $('code', btn.parentNode);
      if (!code || !navigator.clipboard) return;
      navigator.clipboard.writeText(code.textContent).then(function () {
        var was = btn.textContent;
        btn.textContent = 'Copied';
        btn.classList.add('ok');
        setTimeout(function () { btn.textContent = was; btn.classList.remove('ok'); }, 1600);
      });
    });
  });

  /* ----------------------------------------------------------- b-roll clips */
  // The decorative clips carry no `autoplay` attribute: they start when they
  // scroll into view and pause again on the way out, so five of them never sit
  // off-screen decoding at once. With reduced motion, or without an
  // IntersectionObserver, nothing is ever started and each <video> shows its
  // poster instead -- which is exactly the intended fallback, so there is no
  // separate no-JS path to maintain.
  var brolls = $$('[data-broll]');
  if (brolls.length && !reduced && 'IntersectionObserver' in window) {
    var bio = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        var v = e.target;
        if (!e.isIntersecting) { if (!v.paused) v.pause(); return; }
        // preload="none" means nothing is buffered until we ask, so load()
        // first -- otherwise the first play() shows a visible stall.
        if (!v.dataset.loaded) {
          // Not every clip in the set has been shot yet. A missing file would
          // otherwise leave an empty bordered box on the page, so a clip that
          // cannot load takes its container down with it and the layout closes
          // up as though it had never been there.
          //
          // Two mechanisms, because neither is sufficient alone. A failing
          // <source> fires `error` at the <source> itself and `error` does not
          // bubble -- so the listener has to run in the *capture* phase to see
          // it from up here. And when the browser exhausts every <source> it
          // sets NETWORK_NO_SOURCE without firing anything at the <video> at
          // all, which only a check after the fact catches.
          var hide = function () {
            var box = v.closest('.broll');
            if (box) box.hidden = true;
          };
          v.addEventListener('error', hide, { capture: true, once: true });
          v.load();
          setTimeout(function () {
            if (v.networkState === 3 /* NETWORK_NO_SOURCE */) hide();
          }, 1200);
          v.dataset.loaded = '1';
        }
        var playing = v.play();
        // Autoplay can still be refused (low power mode, a strict setting).
        // That leaves the poster up, which is a perfectly good outcome.
        if (playing && playing.catch) playing.catch(function () {});
      });
    }, { threshold: 0.25 });
    brolls.forEach(function (v) { bio.observe(v); });
  }

  // The hero clip is the one that autoplays from the markup, so reduced motion
  // has to stop it here too: `display:none` hides a video without reliably
  // stopping it decoding.
  var heroVideo = $('.hero-video');
  if (heroVideo && reduced) {
    heroVideo.removeAttribute('autoplay');
    heroVideo.pause();
  }

  /* ------------------------------------------------------- screenshot tilt */
  // A few degrees of parallax on the hero shot as the pointer crosses the
  // viewport. Skipped entirely on touch and under reduced motion.
  var tilt = $('.shot-tilt');
  if (tilt && !reduced && matchMedia('(hover: hover)').matches) {
    addEventListener('pointermove', function (e) {
      var x = (e.clientX / innerWidth - 0.5) * 2;
      var y = (e.clientY / innerHeight - 0.5) * 2;
      tilt.style.setProperty('--tilt-y', (x * 2.2).toFixed(2) + 'deg');
      tilt.style.setProperty('--tilt-x', (-y * 1.4).toFixed(2) + 'deg');
    }, { passive: true });
  }
})();
