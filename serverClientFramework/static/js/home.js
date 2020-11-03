$(document).ready(function () {
    ajaxSetup();
    $(".sidenav").sidenav();
});

function ajaxSetup () {
    console.log("Setting up ajax");
    $.ajaxSetup({
        beforeSend: (xhr, settings) => {
            if (!/^(GET|HEAD|OPTIONS|TRACE)$/i.test(settings.type) && !this.crossDomain) {
                xhr.setRequestHeader("X-CSRFToken", $("#csrf-token").val());
            }
        }
    })
}

$("#process-experimental-data").click(function(e) {
    var formData = new FormData();
    formData.append("video", $("#video")[0].files[0]["name"]);
    formData.append("rois", $("#rois-label")[0].files[0]["name"]);
    formData.append("ib4", $("#ib4-label")[0].files[0]["name"]);
    formData.append("cgrp", $("#cgrp-label")[0].files[0]["name"]);
    formData.append("path", $("#folder-path").val());
    $.ajax({
        type: "POST",
        url: document.location.protocol+"//"+document.location.hostname+":"+document.location.port+"/processExperimentalData",
        data: formData,
        cache: false,
        processData: false,
        contentType: false,
        success: (response) => {
            console.log("post success!");
        },
        error: (xhr, errmsg, err) => {
            console.log(xhr, errmsg, err);
        }
    });
});