// Auto-scroll the run log to bottom when the modal is open
(function(){
  function tick(){
    try {
      var modal = document.getElementById('runlog-modal');
      var el = document.getElementById('runlog-content');
      if (!modal || !el) return;
      // Bootstrap modal adds 'show' when visible
      if (modal.classList && modal.classList.contains('show')) {
        el.scrollTop = el.scrollHeight;
      }
    } catch (e) {}
  }
  setInterval(tick, 800);
})();

