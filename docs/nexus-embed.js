/* Runs on bioelectricity.net, next to the knowledge-graph iframe.
   It reads the Nexus login already in this browser and tells the graph who is signed in.
   Paste beside the iframe:
   <iframe id="ee123-graph" data-ee123 src="https://fractastical.github.io/bioelectricity-knowledge-graph/?embed=1" style="width:100%;height:920px;border:0"></iframe>
   <script src="https://fractastical.github.io/bioelectricity-knowledge-graph/nexus-embed.js"></script>
*/
(function () {
  var FRAME_ORIGIN = "https://fractastical.github.io";
  var iframe = document.getElementById("ee123-graph") || document.querySelector("iframe[data-ee123]");
  if (!iframe) return;
  var user = null;

  function storageKey(id) {
    return "ee123_lp_" + id;
  }

  function readPile(id) {
    try {
      var parsed = JSON.parse(localStorage.getItem(storageKey(id)) || "[]");
      return Array.isArray(parsed) ? parsed : [];
    } catch (e) {
      return [];
    }
  }

  async function currentUser() {
    var token = localStorage.getItem("bionexus_auth_token");
    if (!token) return null;
    var res = await fetch("/api/auth/me", { headers: { Authorization: "Bearer " + token } });
    if (!res.ok) return null;
    return res.json();
  }

  async function sendAuth() {
    if (!iframe.contentWindow) return;
    user = await currentUser();
    var verified = user && user.id != null ? readPile(user.id) : [];
    iframe.contentWindow.postMessage({
      type: "BIONEXUS_AUTH",
      user: user ? { id: user.id, username: user.username || "", email: user.email || "" } : null,
      verified: verified
    }, FRAME_ORIGIN);
  }

  window.addEventListener("message", function (event) {
    if (event.origin !== FRAME_ORIGIN) return;
    var data = event.data || {};
    if (data.type === "EE123_READY") sendAuth();
    if (data.type === "EE123_SIGNIN") {
      var next = location.pathname + location.search;
      location.href = "/auth?redirect=" + encodeURIComponent(next);
    }
    if (data.type === "EE123_PROGRESS" && user && String(data.userId) === String(user.id) && Array.isArray(data.verified)) {
      localStorage.setItem(storageKey(user.id), JSON.stringify(data.verified));
    }
  });

  iframe.addEventListener("load", sendAuth);
})();
